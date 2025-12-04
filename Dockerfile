FROM continuumio/miniconda3:latest

# Create conda environment with Python 3.10
RUN conda create -y -n mol python=3.10

# Install RDKit inside the environment
RUN /bin/bash -c "source activate mol && conda install -y -c conda-forge rdkit=2022.03.5"

# Install Python packages
RUN /bin/bash -c "source activate mol && \
    pip install streamlit requests py3Dmol spacy && \
    python -m spacy download en_core_web_sm"

WORKDIR /app
COPY . /app

EXPOSE 8501

CMD ["/bin/bash", "-lc", "source activate mol && streamlit run app.py --server.port=8501 --server.address=0.0.0.0"]
