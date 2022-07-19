import pandas as pd


def parse_known_pathogenic_strs_tsv(known_pathogenic_strs_tsv):
    """Parses a TSV of known pathogenic STR loci into a json data structure.

    Args:
        known_pathogenic_strs_tsv (str): path of the known pathogenic STRs .tsv file

    Returns:
        dict: Returns a dictionary of dictionaries which maps each known pathogenic STR locus id to a dictionary of
         info about that locus - for example, the "FMR1" entry looks like:

            {
            ...
             'FMR1': {
                'Diseases': [
                    {'Inheritance': 'XD',
                        'IntermediateRange': '41-54',
                        'Name': 'Fragile X tremor/ataxia syndrome',
                        'NormalMax': 40,
                        'OMIM': '300623',
                        'PathogenicMin': 55,
                        'Symbol': 'FXTAS'},
                    {'Inheritance': 'XD',
                        'IntermediateRange': '41-54',
                        'Name': 'Fragile X Syndrome',
                        'NormalMax': 40,
                        'OMIM': '300624',
                        'PathogenicMin': 201,
                        'Symbol': 'FXS'}],
                'GeneRegion': "5'-UTR",
                'ReferenceRepeatUnit': 'CGG'},
            ...
            }
    """
    print(f"Loading {known_pathogenic_strs_tsv}")
    known_pathogenic_strs_df = pd.read_table(known_pathogenic_strs_tsv, keep_default_na=False)

    known_pathogenic_strs_df.rename(columns={
        "TR Location (hg38, 1-based)       NOTE: IGV is also 1-based": "TR Location (hg38, 1-based)",
        "TR Location (hg37, 1-based)       NOTE: IGV is also 1-based": "TR Location (hg37, 1-based)",
    }, inplace=True)
    known_pathogenic_strs_df["TR Location (hg38, 1-based)"] = known_pathogenic_strs_df[
        "TR Location (hg38, 1-based)"].str.replace("chr", "")
    known_pathogenic_strs_df["TR Location (hg37, 1-based)"] = known_pathogenic_strs_df[
        "TR Location (hg37, 1-based)"].str.replace("chr", "")

    known_pathogenic_strs_lookup = {}
    for _, row in known_pathogenic_strs_df.iterrows():

        known_pathogenic_strs_lookup[row["LocusId"]] = {
            "Diseases": [],
            "GeneName": row["Gene"],
            "GeneRegion": row["Gene Region"],
            "ReferenceRepeatUnit": row["RU"],
        }

        disease_symbols = row["DiseaseSymbol"].split(";")
        disease_names = row["Disease"].split(";")
        disease_notes = row["DiseaseNotes"].split(";")
        normal_max_thresholds = row["NormalMax gnomAD"].split(";")
        pathogenic_min_thresholds = row["PathogenicMin gnomAD"].split(";")
        omim_disease_links = row["OMIM Disease Link"].split(";")

        field_content_lengths = list(map(len, [
            disease_symbols, disease_names, disease_notes, normal_max_thresholds, pathogenic_min_thresholds,
            omim_disease_links,
        ]))

        if len(set(field_content_lengths)) > 1:
            raise ValueError(f"Different numbers of ';' separated entries found for {row['LocusId']}: {field_content_lengths}")

        intermediate_ranges = [row["IntermediateRange (STRipy)"]] * len(disease_symbols)
        inheritance_modes = [row["InheritanceMode"]] * len(disease_symbols)

        for (
            disease_symbol,
            disease_name,
            disease_note,
            normal_max,
            pathogenic_min,
            intermediate_range,
            inheritance_mode,
            omim_disease_link,
        ) in zip(
            disease_symbols,
            disease_names,
            disease_notes,
            normal_max_thresholds,
            pathogenic_min_thresholds,
            intermediate_ranges,
            inheritance_modes,
            omim_disease_links,
        ):

            disease_info = {
                "Symbol": disease_symbol.strip(),
                "Name": disease_name.strip(),
                "Inheritance": inheritance_mode.strip(),
                "OMIM": omim_disease_link.strip().split("/")[-1],  # OMIM disease id
            }
            disease_note = disease_note.strip()
            normal_max = normal_max.strip()
            pathogenic_min = pathogenic_min.strip()
            intermediate_range = intermediate_range.strip()
            if disease_note:
                disease_info["Note"] = disease_note
            if normal_max:
                disease_info["NormalMax"] = int(normal_max)
            if pathogenic_min:
                disease_info["PathogenicMin"] = int(pathogenic_min)
            if intermediate_range:
                disease_info["IntermediateRange"] = intermediate_range

            known_pathogenic_strs_lookup[row["LocusId"]]["Diseases"].append(disease_info)

    return known_pathogenic_strs_lookup
