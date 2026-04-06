/*
 * Consolidate per-library read statistics across all pipeline steps.
 * Produces a single TSV summary showing reads retained at each step.
 *
 * Works for both v3 (demux + extract + trim + sync) and v2 (trim + sync only).
 */
process consolidate_stats {
    tag "stats ${batch_id}"
    label "pystats"

    publishDir "${params.output_dir}/pipeline_stats", mode: 'copy'

    input:
        val batch_id
        path stat_files   // collected *.tsv from all available steps

    output:
        path "${batch_id}_read_stats.tsv"

    script:
    """
    #!/usr/bin/env python3
    import glob, csv, os

    # ── Collect all TSV fragments ─────────────────────────────────────────────
    rows = []
    for f in glob.glob('*.tsv'):
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter='\\t')
            for row in reader:
                row['_file'] = os.path.basename(f)
                rows.append(row)

    # ── Organise counts by library ────────────────────────────────────────────
    #
    # demux_counts:   run, library, step(raw_r3|demuxed), reads       [v3 only]
    # extract_stats:  library, run, step(demux_ids|extracted), reads   [v3 only]
    # trim_counts:    sample, step(before_trim|after_trim), reads
    # sync_stats:     sample, step(before_trim_sync|after_trim_sync), reads
    #
    # Sample IDs use the pattern: library__run (Kotov) or library (Briggs / v2)

    from collections import defaultdict
    lib_stats = defaultdict(dict)

    for r in rows:
        step = r.get('step', '')
        reads = int(r.get('reads', 0))

        if step == 'raw_r3':
            run = r['run']
            lib_stats['__pool__'][f'raw_r3_{run}'] = reads

        elif step == 'demuxed' and 'run' in r and 'library' in r:
            # from demux_counts only (extract_stats now uses 'demux_ids')
            lib = r['library']
            run = r['run']
            lib_stats[lib][f'demuxed_{run}'] = reads

        elif step == 'demux_ids':
            # from extract_stats — redundant with demuxed, kept for validation
            pass

        elif step == 'extracted':
            lib = r['library']
            run = r['run']
            lib_stats[lib][f'extracted_{run}'] = reads

        elif step in ('before_trim', 'after_trim'):
            sample = r['sample']
            if '__' in sample:
                lib, run = sample.split('__', 1)
                lib_stats[lib][f'{step}_kotov_{run}'] = reads
            else:
                lib_stats[sample][f'{step}_briggs'] = reads

        elif step in ('before_trim_sync', 'after_trim_sync'):
            sample = r['sample']
            if '__' in sample:
                lib, run = sample.split('__', 1)
                lib_stats[lib][f'{step}_kotov_{run}'] = reads
            else:
                lib_stats[sample][f'{step}_briggs'] = reads

    # ── Write summary ─────────────────────────────────────────────────────────
    pool_stats = lib_stats.pop('__pool__', {})

    with open('${batch_id}_read_stats.tsv', 'w') as out:
        out.write('library\\tsource\\tstep\\treads\\tpct_of_previous\\n')

        # Pool-level raw counts
        for key in sorted(pool_stats):
            out.write(f'ALL\\tpool\\t{key}\\t{pool_stats[key]}\\t\\n')

        for lib in sorted(lib_stats):
            s = lib_stats[lib]

            for run in ('runA', 'runB'):
                demuxed = s.get(f'demuxed_{run}')
                extracted = s.get(f'extracted_{run}')
                before_trim = s.get(f'before_trim_kotov_{run}')
                after_trim = s.get(f'after_trim_kotov_{run}')
                after_sync = s.get(f'after_trim_sync_kotov_{run}')

                if demuxed is not None:
                    pool_raw = pool_stats.get(f'raw_r3_{run}', 0)
                    pct = f'{100*demuxed/pool_raw:.1f}%' if pool_raw else ''
                    out.write(f'{lib}\\tkotov_{run}\\tdemuxed\\t{demuxed}\\t{pct}\\n')
                if extracted is not None:
                    pct = f'{100*extracted/demuxed:.1f}%' if demuxed else ''
                    out.write(f'{lib}\\tkotov_{run}\\textracted\\t{extracted}\\t{pct}\\n')
                if before_trim is not None:
                    out.write(f'{lib}\\tkotov_{run}\\tbefore_trim\\t{before_trim}\\t\\n')
                if after_trim is not None:
                    pct = f'{100*after_trim/before_trim:.1f}%' if before_trim else ''
                    out.write(f'{lib}\\tkotov_{run}\\tafter_trim\\t{after_trim}\\t{pct}\\n')
                if after_sync is not None:
                    ref = after_trim or before_trim or 0
                    pct = f'{100*after_sync/ref:.1f}%' if ref else ''
                    out.write(f'{lib}\\tkotov_{run}\\tafter_sync\\t{after_sync}\\t{pct}\\n')

            # Briggs
            before_trim_b = s.get('before_trim_briggs')
            after_trim_b = s.get('after_trim_briggs')
            after_sync_b = s.get('after_trim_sync_briggs')

            if before_trim_b is not None:
                out.write(f'{lib}\\tbriggs\\traw\\t{before_trim_b}\\t\\n')
            if after_trim_b is not None:
                pct = f'{100*after_trim_b/before_trim_b:.1f}%' if before_trim_b else ''
                out.write(f'{lib}\\tbriggs\\tafter_trim\\t{after_trim_b}\\t{pct}\\n')
            if after_sync_b is not None:
                ref = after_trim_b or before_trim_b or 0
                pct = f'{100*after_sync_b/ref:.1f}%' if ref else ''
                out.write(f'{lib}\\tbriggs\\tafter_sync\\t{after_sync_b}\\t{pct}\\n')

            # Total to aligner
            total = 0
            for run in ('runA', 'runB'):
                v = s.get(f'after_trim_sync_kotov_{run}') or s.get(f'after_trim_kotov_{run}', 0)
                if v: total += v
            briggs_final = after_sync_b or after_trim_b or 0
            if briggs_final: total += briggs_final
            if total:
                out.write(f'{lib}\\tcombined\\ttotal_to_aligner\\t{total}\\t\\n')

    print('Stats written to ${batch_id}_read_stats.tsv')
    """
}
