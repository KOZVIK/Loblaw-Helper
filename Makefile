PYTHON ?= python

.PHONY: setup pipeline dashboard test

setup:
	$(PYTHON) -m pip install --requirement requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) analyze_frequencies.py
	$(PYTHON) analyze_response.py
	$(PYTHON) analyze_subsets.py

dashboard:
	$(PYTHON) -m streamlit run dashboard.py

test:
	$(PYTHON) -m unittest discover -s tests -v
