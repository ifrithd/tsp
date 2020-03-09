from tsp import tsp
from PIL import Image, ImageTk
import tkinter as tk
from tkinter import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from pylab import mpl


class City:
    def __init__(self, name, x, y):
        self.name = name
        self.x = float(x)
        self.y = float(y)

class Window:
    def __init__(self):
        self.cities = []
        self.read_cities()

        # 使得可以显示中文
        mpl.rcParams['font.sans-serif'] = ['SimHei']

        self.root = Tk()
        self.root.title('TSP')
        self.root.resizable(width=False, height=False)

        # 画布控件 位于最上方
        figure = Figure(figsize=(6.5, 3.5))
        self.f = figure.add_subplot(1, 1, 1)
        self.canvas = FigureCanvasTkAgg(figure, self.root)
        self.canvas.get_tk_widget().pack()
        self.clean_canvas()

        tk.Label(self.root, text='请选择要经过的城市：').pack()

        # 复选框部分
        self.v = []
        r, c = 0, 0
        sep = Frame()
        sep.pack()
        for i in range(len(self.cities)):
            self.v.append(tk.IntVar())
            b = tk.Checkbutton(sep, text=str(i) + self.cities[i].name, variable=self.v[-1])
            b.grid(row=r, column=c)
            if (i + 1) % 6 == 0:
                r, c = r + 1, 0
            else:
                r, c = r, c + 1
        self.v.append(tk.IntVar())
        b = tk.Checkbutton(sep, text='全选', variable=self.v[-1])
        b.grid(row=r, column=c)

        # 输入框部分
        sep = Frame()
        sep.pack()
        tk.Label(sep, text='请输入初始城市编号：').pack()
        self.textStr = tk.StringVar()
        textEntry = tk.Entry(sep, textvariable=self.textStr)
        self.textStr.set("")
        textEntry.pack()

        # 运行按钮
        tk.Button(sep, text="Run!", command=self._run).pack()

        self.show_dis = tk.Label(sep, text="")
        self.show_dis.pack()
        mainloop()

    def read_cities(self):
        with open('data.txt') as f:
            for line in f:
                self.cities.append(City(*line.split(',')))

    def clean_canvas(self, start_cities=None):
        self.f.clear()
        for i in range(len(self.cities)):
            if start_cities is not None and i == start_cities:
                self.f.scatter(self.cities[i].x, self.cities[i].y, color='b', marker='*')
            else:
                self.f.scatter(self.cities[i].x, self.cities[i].y, color='g', marker='.', alpha=0.4)
            self.f.text(self.cities[i].x, self.cities[i].y, self.cities[i].name, fontsize=7)
        self.canvas.draw()

    def _run(self):
        # 读出选择的城市与起始城市
        selected = [item.get() for item in self.v]
        if selected[-1]:
            selected = [1] * (len(selected)-1)
        else:
            del selected[-1]
        try:
            start_city = int(self.textStr.get())
        except ValueError:
            start_city = -1

        if sum(selected) <= 2:
            self.show_dis['text'] = 'Error： 城市数过少'
            return

        input_cities = []
        for i in range(len(self.cities)):
            if selected[i]:
                input_cities.append(self.cities[i])
                if start_city == -1:
                    start_city = i

        # 计算并绘制图像
        self.clean_canvas(start_city)
        self.show_dis['text'] = '计算中'
        self.root.update()
        ans = tsp(input_cities)
        self.draw_ans(input_cities, ans.dna)
        self.show_dis['text'] = 'distance: {}'.format(ans.distance)

    def draw_ans(self, cities, ans):
        for i in range(len(ans)):
            cityA = cities[ans[i]]
            cityB = cities[ans[(i + 1) % len(ans)]]
            self.f.plot([cityA.x, cityB.x], [cityA.y, cityB.y], color='r')
        self.canvas.draw()


if __name__ == '__main__':
    Window()