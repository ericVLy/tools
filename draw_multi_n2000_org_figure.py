import os
import sys
from tkinter import filedialog
import matplotlib.pyplot as plt


from lib.n2000_txt_draw import DrawFromN2000
# from lib.n2000_txt_to_csv import n2000_txt_to_csv


def get_file_path():
    return filedialog.askdirectory()


if __name__ == "__main__":
    # file_names = os.listdir(get_file_path())
    for filepath,dirnames,filenames in os.walk(get_file_path()):
        for filename in filenames:
            if filename.endswith(".org"):
                file_full_name = os.path.join(filepath, filename)
                # print(file_full_name)
                DrawFromN2000(org_file=file_full_name).draw_from_org_file()
    # plt.show()