"""Read-only support for historical PanTrans GTF inputs."""

from .align_filter import transcript_to_gene_id


def _parse_gtf_attributes(attrs):
    parsed = {}
    for item in attrs.strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.replace('"', "").split(" ", 1)
        parsed[key] = value
    return parsed


def load_gtf_transcript_models(gtf_path):
    """Load validated transcript exon and splice models from an unrenamed GTF."""
    transcript_gene = {}
    transcript_context = {}
    exon_coords = {}
    transcript_features = set()

    with open(gtf_path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.rstrip("\n")
            if not stripped or stripped.startswith("#"):
                continue

            fields = stripped.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{gtf_path}:{line_number}: expected 9 tab-separated GTF fields"
                )

            sequence_name, _, feature, start, end, _, strand, _, attr_text = fields
            if feature not in ("transcript", "exon"):
                continue

            attrs = _parse_gtf_attributes(attr_text)
            missing_attrs = [
                key for key in ("gene_id", "transcript_id") if not attrs.get(key)
            ]
            if missing_attrs:
                raise ValueError(
                    f"{gtf_path}:{line_number}: {feature} is missing "
                    + ", ".join(missing_attrs)
                )

            gene_id = attrs["gene_id"]
            transcript_id = attrs["transcript_id"]
            if sequence_name != gene_id:
                raise ValueError(
                    f"{gtf_path}:{line_number}: sequence name {sequence_name!r} "
                    f"does not match gene_id {gene_id!r}"
                )
            derived_gene_id = transcript_to_gene_id(transcript_id)
            if derived_gene_id != gene_id:
                raise ValueError(
                    f"{gtf_path}:{line_number}: transcript_id {transcript_id!r} "
                    f"maps to gene {derived_gene_id!r}, but gene_id is {gene_id!r}"
                )

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError as exc:
                raise ValueError(
                    f"{gtf_path}:{line_number}: invalid {feature} coordinates "
                    f"{start!r}-{end!r}"
                ) from exc
            if start_i < 1 or end_i < start_i:
                raise ValueError(
                    f"{gtf_path}:{line_number}: invalid {feature} coordinates "
                    f"{start_i}-{end_i}"
                )

            context = (sequence_name, strand)
            previous_context = transcript_context.setdefault(transcript_id, context)
            previous_gene = transcript_gene.setdefault(transcript_id, gene_id)
            if previous_context != context or previous_gene != gene_id:
                raise ValueError(
                    f"{gtf_path}:{line_number}: transcript {transcript_id!r} has "
                    "mixed sequence names or strands"
                )

            if feature == "exon":
                exon_coords.setdefault(transcript_id, []).append((start_i, end_i))
            else:
                transcript_features.add(transcript_id)

    missing_transcript_features = sorted(set(transcript_gene) - transcript_features)
    if missing_transcript_features:
        raise ValueError(
            f"{gtf_path}: transcript {missing_transcript_features[0]!r} has no "
            "transcript feature"
        )

    splice_sites = {}
    sorted_exon_coords = {}
    for transcript_id in transcript_gene:
        coords = sorted(exon_coords.get(transcript_id, []))
        if not coords:
            raise ValueError(
                f"{gtf_path}: transcript {transcript_id!r} has no exon records"
            )

        transcript_splice_sites = []
        for previous_exon, next_exon in zip(coords, coords[1:]):
            if next_exon[0] <= previous_exon[1] + 1:
                raise ValueError(
                    f"{gtf_path}: transcript {transcript_id!r} has overlapping "
                    "or adjacent exons"
                )
            transcript_splice_sites.append(
                (previous_exon[1] + 1, next_exon[0] - 1)
            )

        sorted_exon_coords[transcript_id] = coords
        splice_sites[transcript_id] = transcript_splice_sites

    return {
        "transcript_context": transcript_context,
        "transcript_gene": transcript_gene,
        "splice_sites": splice_sites,
        "exon_coords": sorted_exon_coords,
    }
