"""draw and view one prg file"""
import sys
from tkinter import filedialog
# import matplotlib.pyplot as plt


from lib.n2000_txt_draw import DrawFromN2000
# from lib.n2000_txt_to_csv import n2000_txt_to_csv


def get_file_path():
    """
    get file path
    """
    return filedialog.askopenfilename()


if __name__ == "__main__":
    file_name = get_file_path()
    if file_name.endswith(".txt"):
        DrawFromN2000(txt_file=file_name, show_plt=True).draw_from_txt_file()
    elif file_name.endswith(".org"):
        DrawFromN2000(org_file=file_name, show_plt=True).draw_from_org_file()
    else:
        sys.exit(-1)
    # n2000_txt_to_csv(txt_file=txt_file)
