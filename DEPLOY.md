# Free deployment (the dashboard is a static page)

The dashboard (`dashboard/index.html`) is one self-contained file, so any free static host works.
`python run.py` itself is only for local use; the hosted site serves the already-built page.
Before deploying, run `python run.py --skip-pipeline` once so `dashboard/index.html` is current (it is already included).

## Option A – GitHub Pages (recommended, automatic)
1. Create a GitHub repo and push this folder:
       git init && git add . && git commit -m "CodonRescore project"
       git branch -M main
       git remote add origin https://github.com/<you>/<repo>.git
       git push -u origin main
2. Repo > Settings > Pages > Source: **GitHub Actions**.
3. The included workflow (.github/workflows/pages.yml) deploys on every push.
   Site: https://<you>.github.io/<repo>/

## Option B – Netlify (drag and drop, 1 minute)
Go to app.netlify.com/drop and drag the `dashboard` folder in. Done.
(Or connect the repo; netlify.toml already sets publish = "dashboard".)

## Option C – Cloudflare Pages / Vercel
Connect the repo, set the output/publish directory to `dashboard`, no build command.
