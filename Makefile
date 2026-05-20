.PHONY: help run clean clean-all install

help:
	@echo "🧬 Restriction Enzyme Analyzer - Makefile Commands"
	@echo ""
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  make run           - Analyze plasmid.fasta and generate gel.png and analysis_report.html"
	@echo "  make install       - Install dependencies into .venv"
	@echo "  make serve         - Launch the Streamlit web app"
	@echo "  make clean         - Remove generated output files"
	@echo "  make clean-all     - Remove generated output files and Python cache"
	@echo "  make help          - Show this help message"
	@echo ""
	@echo "Quick Start:"
	@echo "  1. Put your plasmid sequence in plasmid.fasta"
	@echo "  2. Run: make run"
	@echo "  3. Open analysis_report.html in a browser"
	@echo "  4. Or run: make serve"

# Default target
.DEFAULT_GOAL := help

.venv:
	@echo "📦 Creating Python virtual environment..."
	python3 -m venv .venv
	@echo "✓ Virtual environment created"

install: .venv
	@echo "📚 Installing dependencies..."
	@. .venv/bin/activate && pip install --upgrade pip
	@. .venv/bin/activate && pip install -r requirements.txt
	@echo "✓ Dependencies installed"

run: .venv plasmid.fasta
	@echo "🧬 Running restriction enzyme analysis..."
	@. .venv/bin/activate && python find_restriction_enzymes.py --fasta-file plasmid.fasta --report analysis_report.html
	@echo ""
	@echo "✓ Analysis complete"
	@echo "📄 Report: analysis_report.html"
	@echo "🖼️  Gel image: gel.png"

serve: .venv
	@echo "🚀 Launching Streamlit app..."
	@. .venv/bin/activate && streamlit run streamlit_app.py

plasmid.fasta:
	@echo "❌ Error: plasmid.fasta not found!"
	@echo "Please add your plasmid sequence to plasmid.fasta and run make run again."
	@exit 1

clean:
	@echo "🧹 Removing generated output files..."
	@rm -f gel.png analysis_report.html
	@echo "✓ Removed gel.png and analysis_report.html"

clean-all: clean
	@echo "🧹 Removing Python cache files..."
	@rm -rf __pycache__ .pytest_cache
	@echo "✓ Removed cache directories"
