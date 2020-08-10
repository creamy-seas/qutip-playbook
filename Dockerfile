FROM python:3.7.4-buster

ENV PYTHONUNBUFFERED 1

WORKDIR /qutip-simulator

# Setup python
RUN pip install --upgrade pip
COPY requirements.txt .
RUN pip install -r requirements.txt

# Jupyter notebook
RUN pip install jupyter_contrib_nbextensions && jupyter contrib nbextension install
RUN jupyter nbextension enable varInspector/main
RUN pip install autopep8
RUN jupyter nbextension enable code_prettify/autopep8
RUN jupyter nbextension enable hinterland/hinterland
RUN jupyter nbextension enable collapsible_headings/main
RUN jupyter nbextension enable highlight_selected_word/main
RUN jupyter nbextension enable splitcell/splitcell
RUN jupyter nbextension enable toc2/main
RUN jupyter nbextension enable select_keymap/main
RUN jupyter nbextension enable toggle_all_line_numbers/main
RUN jt -t grade3 -cursc o -cursw 5 -T -cellw 80%
COPY ./support_files/notebook.json /root/.jupyter/nbconfig/notebook.json

COPY . .
