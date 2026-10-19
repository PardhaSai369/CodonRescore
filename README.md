# CodonRescore – Gene prediction research project + dashboard

## Run it (one command)
    pip install -r requirements.txt
    python run.py
This runs the full analysis (~1 min), builds the dashboard and opens it at http://localhost:8000.
Use `python run.py --skip-pipeline` to reuse the saved results and just launch the dashboard.

## Layout
- run.py                 launcher (pipeline -> dashboard -> local server)
- code/                  codonrescore.py (method), evaluate.py, analysis.py, export_dashboard.py, gen_tables.py, make_slides.py
- data/                  real NCBI records used
- dashboard/             template.html + data.json -> index.html (also openable directly by double-click)
- paper/                 paper_ieee.tex, paper_springer.tex (shared body.tex); compile with pdflatex twice
- slides/, figs/, results/  presentation, figures, result JSON
