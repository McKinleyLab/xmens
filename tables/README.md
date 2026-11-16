# tables: the QC supplementary table

`supplementary_table_qc.ipynb` builds the QC workbook, `Supplementary_Table_QC.xlsx`, for the final objects: every
library from Cell Ranger to the all-cells and fibroblast/decidual objects, with the nuclei each step removed and why.
It reads the run's own files (Cell Ranger outputs, the ingest and QC logs, the per-sample objects) and the final
objects, and recomputes the per-nucleus QC calls from the counts. Run it in JupyterLab with the `rsc_uv` kernel; it
writes the workbook next to itself.

`check_supplementary_table_qc.py` checks every value of the workbook against the files it reports on, with its own
parsers and its own implementation of the QC rules:

    python check_supplementary_table_qc.py Supplementary_Table_QC.xlsx out_folder

The MI table is made by `mi/` (step 5).
