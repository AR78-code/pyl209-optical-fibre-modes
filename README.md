# PYL209 Assignment 5

Open **Fibre_Demo.html** in a browser for the demonstration. Use **Assignment5_Report.html** to read the derivation and results, or **Assignment5.ipynb** for the runnable notebook.

The demo has two parameters: waveguide length L (0–5 mm) and lens ratio D/f (0.10–1.00). It shows the output in (x,y) and (r,φ), radial profiles, and individual LP mode fields. D remains 20 mm while f changes.

The calculation follows Ghatak, Chapter 13, §§13.5–13.6. It assumes a centred Gaussian beam with field radius 2.5 mm at the lens because the assignment does not specify the incident profile. Six radial modes are excited by this symmetric input. Their amplitudes come from overlap integrals, and their fields are added with the propagation factor exp(−iβL).

## Files

- `fibre.py`: mode roots, field normalization, focusing, excitation and propagation.
- `DERIVATION.md`: derivation with textbook equation references.
- `Assignment5.ipynb`: numerical calculation and plots.
- `validate.py`: numerical checks.
- `results/`: mode data and figures.
- `dashboard.template.html`: source template used to build the demo;
-  Open `Fibre_Demo.html` to run it.
-  Under Combined Near Field, waveguide length L and lens ratio D/f can be varied, and the Transverse Pattern, Near-Field Distribution and Radial Intensity Profile can be seen; both intensity and real field.
-  Under Individual Modes, Modal Pattern, Modal Distribution and Radial Field can be seen for various modes; both intensity and signed field.

## Running the Python calculation

With Python 3.11 or newer, from this folder:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python validate.py
.venv/bin/python build_project.py
.venv/bin/python execute_report.py
```

On Windows, use `.venv\Scripts\python` instead. The HTML demo already contains the calculated data and does not require Python or an internet connection. The report uses MathJax for equations, so its HTML version needs internet access to load that library.

`build_project.py` recreates the plots, demo and notebook. `execute_report.py` runs the notebook and exports the report. `package_project.py` creates the submission ZIP. A separate Matplotlib version is available through `demo.py` if a Python GUI backend is installed.
