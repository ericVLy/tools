"""draw all org files to png and csv in a path"""
import os
# import sys
from tkinter import filedialog
import matplotlib.pyplot as plt


from lib.n2000_txt_draw import DrawFromN2000
from lib.log import logprint as log


def get_file_path():
    """ get file path """
    return filedialog.askdirectory()


if __name__ == "__main__":
    filepath = get_file_path()
    # file_names = os.listdir(filepath)
    for filepath, dirnames, file_names in os.walk(filepath):
        for filename in file_names:
            if filename.endswith(".org"):
                file_full_name = os.path.join(filepath, filename)
                log.info("drawing %s", file_full_name)
                DrawFromN2000(org_file=file_full_name, show_plt=False).draw_from_org_file()
                log.info("finished draw")
                plt.close()
