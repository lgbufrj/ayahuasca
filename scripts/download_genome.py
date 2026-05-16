import os
import requests
from data import GENOME_PATH, organisms
import gzip, lzma
import shutil

OUTPUT_PATH = GENOME_PATH + "/{organism}/{phase}/{type}.{format}"

USERNAME = "enteogenos"
PASSWORD = "#plant_chacmar*"

def decompress_file(filepath, compressed_format, final_format):

    if compressed_format == "gz":

        output_path = filepath[:-3]+f".{final_format}"  # remove .gz

        print(f"Decompressing {filepath} -> {output_path}")

        with gzip.open(filepath, "rb") as f_in:
            with open(output_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove(filepath)

        print("Decompression complete.")

    else:
        print("[download_genome] Unsupported compression format.")

def download_file(organism, file_type, genome_data):

    url = genome_data["url"]

    phase = "phased" if genome_data["phased"] else "non_phased"
    compressed = genome_data["compressed"]

    initial_file_format = "gz" if compressed else genome_data["final_format"]

    output_path = OUTPUT_PATH.format(
        organism=organism,
        phase=phase,
        type=file_type,
        format=initial_file_format
    )

    # If final genome file doesnt exist
    if not os.path.exists(output_path.replace(f".{initial_file_format}", f".{genome_data['final_format']}")):
        
        # If compressed file doesnt exist or is not a compressed genome
        if (compressed and not os.path.exists(output_path)) or (not compressed):

            os.makedirs(os.path.dirname(output_path), exist_ok=True)

            print(f"Downloading {file_type} for {organism}")
            print(f"URL: {url}")
            print(f"Output: {output_path}")

            try:

                headers = {
                    "User-Agent": (
                        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                        "AppleWebKit/537.36 (KHTML, like Gecko) "
                        "Chrome/136.0 Safari/537.36"
                    )
                }

                response = requests.get(
                    url,
                    auth=(USERNAME, PASSWORD),
                    headers=headers,
                    stream=True
                )

                response.raise_for_status()

                with open(output_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)

                print("Download complete.")

            except requests.exceptions.RequestException as e:
                print(f"Failed downloading {file_type} for {organism}")
                print(e)
                
            if compressed:
                decompress_file(output_path, initial_file_format, genome_data["final_format"])
        else:
            print(f"Compressed {file_type} for {organism} already exists. Skipping download and decompressing.")
            decompress_file(output_path, initial_file_format, genome_data["final_format"])
            
    else:
        print(f"{file_type} for {organism} already exists. Skipping download.")


if __name__ == "__main__":

    for organism, organism_data in organisms.items():

        for genome_type, genome_data in organism_data["genome_files"].items():

            if genome_data["download"] is False:
                continue

            download_file(
                organism=organism,
                file_type=genome_type,
                genome_data=genome_data
            )