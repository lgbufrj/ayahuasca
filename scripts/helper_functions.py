import re

# Headers of the caupuri/tucunaca predicted proteomes carry the annotations
# needed downstream (transmembrane domain, signal peptide, chromosome).
CAUPURI_DESCRIPTION_PATTERN = re.compile(
    r"([^|]+)\|([^|]+).*?transmembrane domain:(\w+)\|signal peptide:(\w+).*?\bchr([0-9]+[a-zA-Z]?)\b"
)

# UniProt FASTA headers: "sp|P00509|AAT_ECOLI Aminotransferase OS=..."
UNIPROT_DESCRIPTION_PATTERN = re.compile(
    r"^(?:[a-z]{2,3}\|)?(?P<accession>[A-Za-z0-9]+)\|(?P<entry>\S+)"
)

# Identifier of a "db|ACCESSION|entry" style header; plain ids pass through.
ACCESSION_PATTERN = re.compile(r"^(?:[a-z]{2,3}\|)?(?P<accession>[A-Za-z0-9_.-]+)\|")


def sequence_key(identifier):
    """Canonical id shared by BLAST hit definitions and genome FASTA records.

    "sp|P00509|AAT_ECOLI" -> "P00509"; "CaupuriG00000000004.1|Putative ..." is
    unchanged, so caupuri/tucunaca keep their previous behaviour.
    """
    identifier = (identifier or "").strip()
    match = ACCESSION_PATTERN.match(identifier)
    if match:
        return match.group("accession")
    return identifier.split()[0] if identifier.split() else ""


def _yes_no(value):
    return "y" if str(value).strip().lower() == "yes" else "n"


def format_fasta_description(description, description_only=False):
    description = (description or "").strip()

    match = CAUPURI_DESCRIPTION_PATTERN.match(description)
    if match:
        id_ = match.group(1)
        name = match.group(2)
        td = _yes_no(match.group(3))
        sp = _yes_no(match.group(4))
        chr_ = match.group(5)
    else:
        match = UNIPROT_DESCRIPTION_PATTERN.match(description)
        if match:
            id_ = match.group("accession")
            name = match.group("entry")
            td = sp = chr_ = ""
        elif "|" in description:
            id_, _, rest = description.partition("|")
            id_ = id_.strip()
            name = rest.split("|")[0].strip()
            td = sp = chr_ = ""
        else:
            id_ = description.split()[0] if description.split() else ""
            name = td = sp = chr_ = ""

    parts = [id_]
    if name:
        parts.append(name)
    if td or sp:
        parts.append(f"td {td} | sp {sp}")
    if chr_:
        parts.append(chr_)
    new_description = " | ".join(part for part in parts if part)

    if description_only:
        return new_description
    return new_description, id_, name, td, sp, chr_