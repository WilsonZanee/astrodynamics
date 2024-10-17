import argparse
import os
import sys
from pathlib import Path

from reportlab.pdfgen import canvas
from pypdf import PdfMerger
from fpdf import FPDF


def merge_pdfs(files):
    dir = os.path.dirname(files[0])
    basename = Path(files[0]).stem
    try:
        Path(os.path.join(dir, "pdf_output")).mkdir()
    except FileExistsError:
        pass
    pdfs = []
    for file in files:
        ext = Path(file).suffix
        if ext == '.pdf':
            pdfs.append(file)
        elif ext == '.py':
            pdfs.append(convert_to_pdf(file))
        else:
            print(f"Error: This program does not support {ext} files.")
            sys.exit(1)

    merger = PdfMerger()
    for pdf in pdfs:
        merger.append(pdf)
    pdf_path = os.path.join(dir, "pdf_output", f"{basename}_merged.pdf")
    print(pdf_path)
    merger.write(pdf_path)

def convert_to_pdf(file):
    dir = os.path.dirname(file)
    basename = os.path.basename(file)
    pdf_path = os.path.join(dir, "pdf_output", f"{basename}.pdf")
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size=12)
 
    with open(file, 'r') as python_file:
        content = python_file.read()
        pdf.multi_cell(0, 10, content)
 
    pdf.output(pdf_path)
    return pdf_path

"""
def convert_to_pdf(file):
    dir = os.path.dirname(file)
    basename = os.path.basename(file)
    pdf_path = os.path.join(dir, f"{basename}.pdf")
    pdf_canvas = canvas.Canvas(pdf_path)
 
    with open(file, 'r') as python_file:
        content = python_file.read()
        pdf_canvas.drawString(100, 800, content)
 
    pdf_canvas.save()
    return pdf_path
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+')
    
    args = parser.parse_args()
    merge_pdfs(args.files)