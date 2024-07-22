import os
import sys
from tkinter import filedialog
import matplotlib.pyplot as plt

from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import Paragraph, SimpleDocTemplate, Image, Table
from reportlab.platypus import Spacer
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics.charts.legends import Legend
from reportlab.lib import  colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm

from lib.n2000_txt_draw import DrawFromN2000
from lib.log import logprint as log


def get_file_path():
    return filedialog.askdirectory()


if __name__ == "__main__":
    filepath = get_file_path()
    file_names = os.listdir(filepath)
    for filename in file_names:
        if filename.endswith(".org"):
            file_full_name = os.path.join(filepath, filename)
            png_file = filename.split(".org")[0]+".png"
            csv_file = filename.split(".org")[0]+"_org.csv"
            if not (os.path.isfile(os.path.join(filepath,  png_file)) and os.path.isfile(os.path.join(filepath,  csv_file))):
                log.info("drawing %s", file_full_name)
                DrawFromN2000(org_file=file_full_name, show_plt=False).draw_from_org_file()
                log.info("finished draw")
                plt.close()
            # build report