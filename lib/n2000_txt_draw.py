# N2000
import argparse
import os
import json
import csv
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from matplotlib.backend_bases import MouseButton
import numpy as np
from numpy import trapz
from scipy.integrate import simpson

sys.path.append(os.path.abspath(".."))
sys.path.append(os.path.abspath("."))
from lib.log import logprint as log



class DrawFromN2000:
    def __init__(self, txt_file=None, org_file=None, export_filename="exp.png", encoding="gbk", show_plt=False):
        self.text_file = txt_file
        self.org_file_name = org_file
        self.encoding = encoding
        self.x_data = None
        self.y_data = None
        self.dy = None
        self.binding_id = None
        self.ax = None
        self.fig = None
        self.line = None
        self.ann = None
        self.event_point = None
        self.heigh_points = []
        self.using_text_file = False
        self.using_org_file = False
        self.show_plt = show_plt

    def __on_move(self, event):
        if event.inaxes == self.ax:
            # log.info(f'data coords {event.xdata} {event.ydata},',
            #     f'pixel coords {event.x} {event.y}')
            x, y = event.xdata, event.ydata
            data = np.interp(x, self.x_data, self.y_data)
            if self.ann is not None:
                self.ann.remove()
            self.ann = self.ax.annotate(f'({x:.2f}, {data:.2f})',xy=(x,y))
            if self.event_point is not None:
                self.event_point.remove()
            self.event_point, = self.ax.plot(x, data, marker='o', color="blue", markersize=5)
            self.fig.canvas.draw_idle()


    def __on_right_click(self, event):
        if event.button is MouseButton.RIGHT:
            log.info('connecting callback')
            if self.binding_id is not None:
                plt.disconnect(self.binding_id)
            self.binding_id = self.fig.canvas.mpl_connect('motion_notify_event', self.__on_move)

    def __on_click(self, event):
        if event.button is MouseButton.LEFT:
            log.info('disconnecting callback')
            plt.disconnect(self.binding_id)
            self.binding_id = None

    def __trapz_area(self, start_index:int, end_index:int):
        area = trapz(y=self.y_data[start_index: end_index+1], dx=100)
        return round(area, 4)

    def __simpson_area(self, start_index:int, end_index:int):
        area = simpson(self.y_data[start_index: end_index+1], dx=100)
        return round(area, 4)

    def __get_heigh_point(self, max_minutes=27):
        self.heigh_points = []
        default_l = 0
        default_h = 0
        default_y = 0
        count = 0
        ndight = 1
        trapz_area_list = []
        simpson_area_list = []
        down = False
        up = False
        for i, y_idx in enumerate(self.y_data):
            if self.x_data[i] > max_minutes:
                break
            if i> 0 and round(y_idx, ndight) > round(self.y_data[i-1], ndight):
                up = True
            elif i> 0 and round(y_idx, ndight) < round(self.y_data[i-1], ndight):
                down = True
                up = False
                continue
            if down and up:
                default_h = i
                if default_h > default_l:
                    y_heigh = np.max(self.y_data[default_l: default_h+1])
                    x_heigh = [x for j,x in enumerate(self.x_data[default_l: default_h+1]) if self.y_data[default_l: default_h+1][j] == y_heigh][0]
                    count+=1
                    trapz_area_list.append(self.__trapz_area(default_l, default_h))
                    simpson_area_list.append(self.__simpson_area(default_l, default_h))
                    self.heigh_points.append({
                                              "峰号": count,
                                              "峰名": None,
                                              "保留时间": x_heigh, 
                                            #   "宽度": round(self.x_data[default_h] - self.x_data[default_l], 4),
                                              "峰高": y_heigh, 
                                              "trapz 峰面积": self.__trapz_area(default_l, default_h),
                                              "trapz 计算含量": 0,
                                              "simpson 峰面积": self.__simpson_area(default_l, default_h),
                                              "simpson 计算含量": 0,
                                              })
                    self.ax.plot(x_heigh, y_heigh, marker='o', color="red", markersize=3)
                    self.ax.annotate(f'{x_heigh:.2f}',xy=(x_heigh,y_heigh), xytext=(x_heigh,y_heigh*1.05), fontsize=8)
                    self.fig.canvas.draw_idle()
                    default_h = i
                    default_l = i
                    default_y = y_idx
                    down = False
        # export csv
        sum_trapz_area = sum(trapz_area_list)
        sum_simpson_area = sum(simpson_area_list)
        for hp in self.heigh_points:
            hp["trapz 计算含量"] = hp["trapz 峰面积"] / sum_trapz_area
            hp["simpson 计算含量"] = hp["simpson 峰面积"] / sum_simpson_area
        if not self.using_org_file and not self.using_text_file and not os.path.exists("export"):
            os.makedirs("export")
            file_name = "export.csv"
            file_path = os.path.join("export", file_name)
        if self.using_text_file:
            file_name = os.path.split(self.text_file)[-1].split(".")[0]+"_text.csv"
            file_path = os.path.join(os.path.dirname(self.text_file), file_name)
        if self.using_org_file:
            file_name = os.path.split(self.org_file_name)[-1].split(".")[0]+"_org.csv"
            file_path = os.path.join(os.path.dirname(self.org_file_name), file_name)
        with open(os.path.join(file_path), "w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=self.heigh_points[0].keys())
            writer.writeheader()
            writer.writerows(self.heigh_points)
        return

    def __save_x_y_data(self, file_name="data.csv"):
        # np.savetxt(os.path.join("export", file_name), np.array([self.x_data, self.y_data]), delimiter=',')
        with open(file_name, "w", newline="") as file:
            fieldnames = ['time','elec']
            data = []
            for i in range(len(self.x_data)):
                data.append({"time": self.x_data[i],
                            "elec": self.y_data[i]})
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(data)

    def draw_from_org_file(self):
        self.using_org_file = True
        self.using_text_file = False
        if self.org_file_name:
            with open(self.org_file_name, "rb") as org_file:
                org = org_file.read()
                point_count = np.frombuffer(org[:4],dtype=np.uint32)[0]
                data_32 = np.frombuffer(org[4:point_count*4+4], dtype=np.int32)
                self.y_data = np.array([round(float(i/1000), 3) for i in data_32[1:]])
                self.x_data = np.array([round(float(i* 1/600), 4) for i in range(1, point_count)])  # 0.01s == 1/600min
                # log.info(list(self.x_data))
                # log.info(list(self.y_data))
                file_name = os.path.split(self.org_file_name)[-1].split(".")[0]+"_org_data.csv"
                file_path = os.path.join(os.path.dirname(self.org_file_name), file_name)
                self.__save_x_y_data(os.path.join(file_path))
                self.__draw()
            return
        assert self.org_file_name, "org file input is None"

    def draw_from_txt_file(self):
        self.using_org_file = False
        self.using_text_file = True
        with open(self.text_file, "r", encoding=self.encoding) as atxt:
            a = atxt.readlines()
            all_data = []
            x_starts = None
            y_sum = None
            index = 0

            for l in range(2, len(a)):
                exp = [string for string in a[l].split("\n")[0].split(" ") if string]
                x_round = float(exp[0])
                y_idx = float(exp[1])
                all_data.append([x_round, y_idx])
            x_data = np.array([x for x, y in all_data])
            for x, y in all_data:
                assert x >= 0, "time shall greater than 0"
            y_data = np.array([y for x, y in all_data])
            assert len(x_data) == len(y_data)
        self.x_data = x_data
        self.y_data = y_data
        # log.info(list(self.x_data))
        # log.info(list(self.y_data))
        file_name = os.path.split(self.text_file)[-1].split(".")[0]+"_text_data.csv"
        file_path = os.path.join(os.path.dirname(self.text_file), file_name)
        self.__save_x_y_data(os.path.join(file_path))
        self.__draw()

    def __draw(self):

        matplotlib.rcParams['font.sans-serif'] = ['SimHei']
        # Y轴范围
        plt.ylim((0, 320))
        plt.xlim((0, np.max(self.x_data)))
        # x轴
        self.fig, self.ax = plt.subplots(figsize=(8, 4))
        self.line, = self.ax.plot(self.x_data, self.y_data, linewidth=0.5)
        self.binding_id = self.fig.canvas.mpl_connect('motion_notify_event', self.__on_move)
        x_major_locator=MultipleLocator(2)
        self.ax.xaxis.set_major_locator(x_major_locator)

        plt.xlabel("时间/min")
        plt.ylabel("电压/mv")
        self.__get_heigh_point()
        if self.using_text_file:
            file_name = os.path.split(self.text_file)[-1].split(".")[0]+".png"
            file_path = os.path.join(os.path.dirname(self.text_file), file_name)
        if self.using_org_file:
            file_name = os.path.split(self.org_file_name)[-1].split(".")[0]+".png"
            file_path = os.path.join(os.path.dirname(self.org_file_name), file_name)
        elif not os.path.exists("export") and not self.using_text_file:
            os.makedirs("export")
            file_path = os.path.join("export", f"{os.path.split(self.text_file)[-1].split(".")[0]}.png")
        self.ax.set_title(file_name.split(".")[0])
        
        plt.savefig(file_path, dpi=300)
        # self.binding_id = plt.connect('motion_notify_event', self.__on_move)
        self.fig.canvas.mpl_connect("button_press_event", self.__on_right_click)
        self.fig.canvas.mpl_connect('button_press_event', self.__on_click)
        if self.show_plt:
            plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='draw n2000 export')
    parser.add_argument('--text_file', "-t", type=str, help='file name', default=r"")
    parser.add_argument('--org_file', "-o", type=str, help='file name', default=r"export\2024-07-19-01-04.org")
    parser.add_argument('--export_filename', "-e", type=str, default=f"exp.png", help='export picture name')

    args = parser.parse_args()
    txt_file = args.text_file
    
    # DrawFromN2000(txt_file=txt_file, org_file=args.org_file, export_filename=args.export_filename).draw_from_txt_file()
    draw = DrawFromN2000(txt_file=txt_file, org_file=args.org_file, export_filename=args.export_filename)
    np.set_printoptions(threshold=sys.maxsize)
    draw.draw_from_org_file()
    # draw.draw_from_txt_file()