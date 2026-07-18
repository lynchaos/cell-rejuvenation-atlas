# Rejuvenation Explorer

Interactive Streamlit app over the pipeline outputs — the "data exploration tool" deliverable.

```bash
# after running any pipeline module:
streamlit run app/streamlit_app.py -- --results results/
```

Tabs correspond 1:1 to analysis modules. All plots read only files under `results/` — the app never recomputes science, keeping it fast and safe for collaborators.

## Deploy (optional)

* **Streamlit Community Cloud**: point it at this repo, entry `app/streamlit_app.py`, and ship a small demo `results/` snapshot in `app/cache/`.
* **AWS App Runner / ECS**: `docker/Dockerfile.app` packages the app for a shareable internal URL.
