"""Create a portable submission archive without environments or source PDFs."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

root=Path(__file__).resolve().parent
names=["README.md","DERIVATION.md","requirements.txt","requirements-lock.txt",
       "fibre.py","demo.py","validate.py","build_project.py","execute_report.py",
       "package_project.py","dashboard.template.html","Fibre_Demo.html",
       "Assignment5.ipynb","Assignment5_Report.html"]
files=[root/name for name in names]
result_names=["modal_distributions.png", "length_patterns.png", "lens_patterns.png",
              "modes.csv", "coupling_sweep.csv", "validation.json", "summary.json"]
files.extend(root/"results"/name for name in result_names)
with ZipFile(root/"Assignment5_Project.zip","w",ZIP_DEFLATED) as z:
    for p in files:
        z.write(p,Path("Assignment5_Project")/p.relative_to(root))
with ZipFile(root/"Assignment5_Project.zip") as z:
    assert z.testzip() is None
    print(f"Verified {len(z.namelist())} files in Assignment5_Project.zip")
