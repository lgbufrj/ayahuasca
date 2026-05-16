import re

def format_fasta_description(description, description_only=False):
    match = re.match(
            r"([^|]+)\|([^|]+).*?transmembrane domain:(\w+)\|signal peptide:(\w+).*?\bchr([0-9]+[a-zA-Z]?)\b",
            description,
        )
    if match:
        id_ = match.group(1)
        name = match.group(2)
        td = "y" if match.group(3).lower() == "yes" else "n"
        sp = "y" if match.group(4).lower() == "yes" else "n"
        chr_ = match.group(5)
        
        new_description = f"{id_} | {name} | td {td} | sp {sp} | {chr_}"
        # .replace(" ", "_")

    if description_only:
        return new_description
    return new_description, id_, name, td, sp, chr_