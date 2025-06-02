#!/usr/bin/python
import os
import glob
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_initial_plotting(options):

    output_dir = os.path.join(options.result_files_directory, 'Reports')

    # Plotte die performance matrix mit den MCC Werten
    plotting_performance(options, output_dir)    

    # Plotter für die presence absence matrix histogramme
    plotting_matrix_histogram(options, output_dir)
    



    return

def _prepare_output_dir(report_directory, output_dir):
    if output_dir is None:
        output_dir = os.path.join(report_directory, 'Reports')
    os.makedirs(output_dir, exist_ok=True)
    return output_dir
    


def _run_rscript(r_script, tsv_file, flag):
    """
    Hilfsfunktion, die das R-Skript mit den gegebenen Argumenten aufruft.
    Gibt ein Tupel (tsv_file, flag, returncode) zurück.
    """
    result = subprocess.run(
        ["Rscript", r_script, tsv_file, flag],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    return (tsv_file, flag, result.returncode)

def plotting_matrix_histogram(options, tsv_dir):
    """
    Für jede TSV-Datei unter tsv_dir, deren Name mit "neighborhood_confusion_" beginnt,
    ruf das R-Skript zweimal auf: einmal mit filter=TRUE und einmal mit filter=FALSE.
    Parallelisiert über ThreadPoolExecutor. Prüft vor jedem Task, ob die jeweilige Output-PDF
    bereits existiert und überspringt den Task in diesem Fall.
    """
    # 1) Pfad zum R-Skript
    r_script = os.path.join(options.plotting_Rscripts, "Plot_confusion_matrix_TP.R")

    # 2) Rekursiv nach allen "neighborhood_confusion_*.tsv"-Dateien suchen
    pattern = os.path.join(tsv_dir, "**", "neighborhood_confusion_*.tsv")
    tsv_files = glob.glob(pattern, recursive=True)

    if not tsv_files:
        print(f"No files matching 'neighborhood_confusion_*.tsv' found in {tsv_dir}.")
        return

    # 3) Tasks erzeugen: für jede Datei zweimal (TRUE und FALSE), aber nur wenn Output noch nicht existiert
    tasks = []
    for tsv_file in sorted(tsv_files):
        base = os.path.splitext(tsv_file)[0]
        for flag in ["TRUE", "FALSE"]:
            # Output-PDF-Name: <basename>_<flag>.pdf
            output_pdf = f"{base}_{flag}.pdf"
            if os.path.exists(output_pdf):
                print(f"[SKIP] '{output_pdf}' existiert bereits – überspringe.")
            else:
                tasks.append((tsv_file, flag))

    if not tasks:
        print("Alle Outputs existieren bereits, keine Tasks gestartet.")
        return

    # 4) Parallel ausführen mit ThreadPoolExecutor
    max_workers = min(len(tasks), os.cpu_count() or 1)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {
            executor.submit(_run_rscript, r_script, tsv_file, flag): (tsv_file, flag)
            for tsv_file, flag in tasks
        }

        for future in as_completed(future_to_task):
            tsv_file, flag = future_to_task[future]
            try:
                tsv_file, flag, returncode = future.result()
                if returncode != 0:
                    print(f"[ERROR] R-Skript schlug fehl für '{tsv_file}' (filter={flag}).")
                else:
                    print(f"[SAVE] Plot erstellt für '{tsv_file}' (filter={flag}).")
            except Exception as exc:
                print(f"[ERROR] Beim Ausführen von '{tsv_file}' (filter={flag}) ist ein Fehler: {exc}")


                
def plotting_performance(options, tsv_dir):
    """
    Sucht rekursiv nach jeder 'all_performance.txt' in tsv_dir, 
    ruft für jede den R-Skript 'Plot_performance_matrix.R' auf und 
    erzeugt eine PDF an derselben Stelle.
    """
    # TODO erweitern um andere Messungen abzubilden
    
    # 1) Pfad zum R-Skript (r_script muss genau so namensgeladen werden)
    r_script = os.path.join(options.plotting_Rscripts, "Plot_performance_matrix.R")

    # 2) Rekursiv nach allen "all_performance.txt"-Dateien suchen
    pattern = os.path.join(tsv_dir, "**", "all_performance.txt")
    tsv_files = glob.glob(pattern, recursive=True)

    if not tsv_files:
        print(f"No files matching 'all_performance.txt' found in {tsv_dir}.")
        return

    # 3) Für jede gefundene TSV-Datei einmal Plot aufrufen
    for tsv_file in sorted(tsv_files):
        # PDF-Ausgabepfad: gleiche Location, aber Endung .pdf
        base, _ = os.path.splitext(tsv_file)
        output_pdf = f"{base}.pdf"

        print(f"Processing performance plot for '{tsv_file}' …")
        result = subprocess.run(
            ["Rscript", r_script, tsv_file, output_pdf, "MCC"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        if result.returncode != 0:
            print(f"[ERROR] R script failed on '{tsv_file}'.")
        else:
            print(f"[OK] PDF created: {output_pdf}")


