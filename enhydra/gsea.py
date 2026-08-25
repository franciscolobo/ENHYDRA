def run_gsea(
    anchor2mean_path: str,
    results_dir: str,
    gene_sets: str,
    organism: str | None = None,
    sources: list[str] | None = None,
    permutations: int = 1000,
    min_size: int = 5,
    max_size: int = 500,
    seed: int = 42,
    fdr_threshold: float = 0.25,
    weight: float = 0,
) -> gp.Prerank:
    """Run GSEApy prerank, optionally fetching annotations from g:Profiler."""
    os.makedirs(results_dir, exist_ok=True)
    ranked = _load_ranked(anchor2mean_path)

    if organism is not None:
        gmt_path = os.path.join(results_dir, "gprofiler_%s.gmt" % organism)

        if os.path.isfile(gmt_path):
            logger.info("GMT file already exists, skipping API call: %s", gmt_path)
        else:
            build_gmt_from_gprofiler(
                gene_ids=ranked["gene_id"].tolist(),
                organism=organism,
                gmt_path=gmt_path,
                sources=sources,
            )

        gene_sets = gmt_path

    # Preflight: fail early with a clear message if IDs don't match
    _check_overlap(ranked, gene_sets)

    logger.info(
        "Running GSEApy prerank: %d genes, gene sets: %s, permutations: %d",
        len(ranked),
        gene_sets,
        permutations,
    )

    results = gp.prerank(
        rnk=ranked,
        gene_sets=gene_sets,
        outdir=results_dir,
        permutation_num=permutations,
        min_size=min_size,
        max_size=max_size,
        seed=seed,
        weight=weight,
        graph_num=0,
        verbose=False,
    )

    sig_rows = results.res2d[
        results.res2d["FDR q-val"] < fdr_threshold
    ]
    n_sig = len(sig_rows)

    logger.info(
        "GSEA complete. %d significant gene sets (FDR < %s).",
        n_sig,
        fdr_threshold,
    )

    if n_sig > 0:
        logger.info(
            "Generating enrichment plots for %d significant gene sets...",
            n_sig,
        )

        plots_dir = os.path.join(results_dir, "prerank")
        os.makedirs(plots_dir, exist_ok=True)

        for term in sig_rows["Term"]:
            safe_term_name = term.replace("/", "_").replace("\\", "_")
            output_path = os.path.join(
                plots_dir,
                f"{safe_term_name}.png",
            )

            gp.gseaplot(
                rank_metric=results.ranking,
                term=term,
                **results.results[term],
                ofname=output_path,
            )

    return results
