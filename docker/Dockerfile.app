FROM python:3.11-slim

WORKDIR /app
COPY app/streamlit_app.py ./streamlit_app.py
COPY app/cache/ ./results/
RUN pip install --no-cache-dir streamlit==1.36.* plotly==5.22.* pandas==2.2.* pyarrow==16.*

EXPOSE 8501
ENTRYPOINT ["streamlit", "run", "streamlit_app.py", "--server.address=0.0.0.0", "--", "--results", "results/"]
