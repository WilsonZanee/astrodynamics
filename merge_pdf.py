import argparse
import os
import sys
from pathlib import Path

from pypdf import PdfMerger
from fpdf import FPDF
import win32com.client
from subprocess import Popen


def merge_pdfs(files):
    dir = os.path.dirname(files[0])
    basename = Path(files[0]).stem
    try:
        Path(os.path.join(dir, "pdf_output")).mkdir()
    except FileExistsError:
        pass
    pdfs = []
    for file in files:
        if file.endswith('.pdf'):
            pdfs.append(os.path.abspath(file))
        elif file.endswith('.py'):
            pdfs.append(py_to_pdf(file))
        elif file.endswith(('.doc', '.docx', 'ppt', 'pptx')):
            pdfs.append(doc_ppt_to_pdf(file))
        else:
            print(f"Cannot convert {file} to pdf - file type unsupported")

    merger = PdfMerger()
    print(pdfs)
    for pdf in pdfs:
        merger.append(pdf)
    pdf_path = os.path.join(dir, "pdf_output", f"merged.pdf")
    print(pdf_path)
    merger.write(pdf_path)

def py_to_pdf(file):
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

def doc_ppt_to_pdf(file):
    file = Path(os.path.abspath(file))
    dir = os.path.dirname(file)
    output_dir = os.path.join(dir, "pdf_output")
    basename = Path(os.path.basename(file)).stem
    pdf_path = os.path.join(dir, "pdf_output", f"{basename}.pdf")
    LIBRE_OFFICE = r"C:\Program Files\LibreOffice\program\soffice.exe"
    p = Popen([LIBRE_OFFICE, '--headless', '--convert-to', 'pdf', '--outdir',
               output_dir, file])
    print([LIBRE_OFFICE, '--convert-to', 'pdf', file])
    p.communicate()
    print(pdf_path)
    return pdf_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+')
    args = parser.parse_args()

    if len(args.files) == 1:
        try:
            files = os.listdir(args.files[0])
            abs_files = []
            for file in files:
                print(os.path.abspath(args.files[0]))
                print(file)
                abs_file = os.path.join(os.path.abspath(args.files[0]), file)
                abs_files.append(abs_file)
                print(abs_files)
            print(abs_files)
            merge_pdfs(abs_files)
        except:
            pass
    else:
        merge_pdfs(args.files)