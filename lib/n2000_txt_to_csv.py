# N2000 txt to csv
import csv
import argparse

import os

def n2000_txt_to_csv(txt_file: str, export_filename="exp.csv"):

    with open(txt_file, "r", encoding="gbk") as atxt:
        a = atxt.readlines()
        all_data = []
        x_starts = None
        y_sum = None

        for l in range(len(a)):
            exp = [string for string in a[l].split("\n")[0].split(" ") if string]
            all_data.append(exp)
    if not os.path.exists("export"):
        os.makedirs("export")
    with open(os.path.join("export", export_filename), "w", encoding="gbk", newline='') as f:
        writer = csv.writer(f)
        writer.writerows(all_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='draw n2000 export')
    parser.add_argument('--file', "-f", type=str, help='file name', required=True)
    parser.add_argument('--export_filename', "-e", type=str, default="exp.csv", help='export file name')
    args = parser.parse_args()
    txt_file = args.file
    n2000_txt_to_csv(txt_file=txt_file, export_filename=args.export_filename)