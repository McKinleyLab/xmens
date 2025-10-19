import argparse
import pandas as pd
from samap.mapping import SAMAP
from samap.utils import save_samap

def run_samap(species1, species2, data1, data2, maps_path, out_path):
    """
    Function to run SAMap integration.

    Parameters:
        species1 (str): Abbreviation for the first species (e.g., 'hs').
        species2 (str): Abbreviation for the second species (e.g., 'mm').
        data1 (str): Path to the data file for the first species.
        data2 (str): Path to the data file for the second species.
        maps_path (str): Path to the BLAST/HMMER output directory.
        out_path (str): Path to save the SAMap object.
    """
    print("Starting SAMap integration...")

    #lLoad in raw data files
    print(f"Species 1 ({species1}) data path: {data1}")  
    print(f"Species 2 ({species2}) data path: {data2}")
    filenames = {species1: data1, species2: data2}

    # create the SAMap object
    print(f"Initializing SAMap with maps path: {maps_path}")
    sm = SAMAP(
        filenames,
        f_maps=maps_path,
        save_processed=True  # Save processed results to `*_pr.h5ad`
    )

    # run SAMap integration
    print("Running SAMap integration with pairwise stitching...")
    sm.run(pairwise=True)

    # access the object with the integrated data
    print("SAMap integration complete. Accessing SAMap object...")
    samap = sm.samap  # SAMap object with integrated data

    # save the SAMap object
    print(f"Saving SAMap object to: {out_path}")
    save_samap(sm, out_path)
    print("SAMap object saved successfully.")
    print("SAMap integration pipeline completed.")

if __name__ == "__main__":
    # set up argument parsing
    parser = argparse.ArgumentParser(description="Run SAMap integration for two species.")
    parser.add_argument("--species1", required=True, help="Abbreviation for the first species (e.g., 'hs').")
    parser.add_argument("--species2", required=True, help="Abbreviation for the second species (e.g., 'mm').")
    parser.add_argument("--data1", required=True, help="Path to the data file for the first species.")
    parser.add_argument("--data2", required=True, help="Path to the data file for the second species.")
    parser.add_argument("--maps", required=True, help="Path to the BLAST/HMMER output directory.")
    parser.add_argument("--out", required=True, help="Path to save the SAMap object.")

    # parse args
    args = parser.parse_args()

    # run SAMap with parsed args
    run_samap(args.species1,
              args.species2,
              args.data1,
              args.data2,
              args.maps,
              args.out)