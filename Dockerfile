FROM continuumio/miniconda

ENV BK_VERSION=2.4.0
ENV PY_VERSION=3.9
ENV NUM_PROCS=2
ENV BOKEH_RESOURCES=cdn
ENV BOKEH_LOG_LEVEL=debug
copy . /app
WORKDIR /app

RUN apt-get install git bash
RUN pip install -r /app/requirements.txt

RUN conda config --append channels bokeh
RUN conda config --append channels districtdatalabs
RUN conda install --yes --quiet python=${PY_VERSION} pyyaml jinja2 bokeh=${BK_VERSION} numpy "nodejs>=8.8" pandas scipy
RUN conda clean -ay

EXPOSE 5006
EXPOSE 80


CMD bokeh serve \
    --allow-websocket-origin="*" \
    --num-procs=${NUM_PROCS} \
    --show /app/portal_hdf5

    