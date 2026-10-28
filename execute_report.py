"""Execute the report with a project-local kernel spec and export HTML."""
import os
import sys
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
kernel=ROOT/"tmp/jupyter/kernels/python3"
kernel.mkdir(parents=True,exist_ok=True)
(kernel/"kernel.json").write_text(json.dumps({"argv":[sys.executable,"-m","ipykernel_launcher","-f","{connection_file}"],"display_name":"Project Python","language":"python"}))
os.environ["JUPYTER_PATH"]=str(ROOT/"tmp/jupyter")
os.environ["JUPYTER_RUNTIME_DIR"]=str(ROOT/"tmp/jupyter/runtime")
os.environ["IPYTHONDIR"]=str(ROOT/"tmp/ipython")
os.environ["MPLCONFIGDIR"]=str(ROOT/".mplconfig")
import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter
nb=nbformat.read(ROOT/"Assignment5.ipynb",as_version=4)
client=NotebookClient(nb,timeout=180,kernel_name="python3",resources={"metadata":{"path":str(ROOT)}})
client.execute()
nbformat.write(nb,ROOT/"Assignment5.ipynb")
html,_=HTMLExporter().from_notebook_node(nb)
html=html.replace("<title>Notebook</title>","<title>Assignment 5 — derivation and results</title>")
html=html.replace("</head>","<style>.jp-RenderedMarkdown > table {margin:1em 0; min-width:500px;} .jp-RenderedMarkdown > table th,.jp-RenderedMarkdown > table td {text-align:left!important; padding:8px 16px!important;} .jp-RenderedMarkdown > table th:first-child,.jp-RenderedMarkdown > table td:first-child {min-width:230px;}</style></head>")
(ROOT/"Assignment5_Report.html").write_text(html)
print("Executed all report cells and exported Assignment5_Report.html")
