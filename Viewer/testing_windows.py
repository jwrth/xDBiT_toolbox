# import tkinter as tk
#
# root = tk.Tk()
#
# canvas1 = tk.Canvas(root, width=400, height=300)
# canvas1.pack()
#
# entry1 = tk.Entry(root)
# canvas1.create_window(200, 140, window=entry1)
#
#
# def getSquareRoot():
#     x1 = entry1.get()
#
#     label1 = tk.Label(root, text=float(x1) ** 0.5)
#     canvas1.create_window(200, 230, window=label1)
#
#
# button1 = tk.Button(text='Get the Square Root', command=getSquareRoot)
# canvas1.create_window(200, 180, window=button1)
#
# root.mainloop()

# import tkinter as tk
#
# top = tk.Tk()
# L1 = tk.Label(top, text="User Name")
# L1.pack( side = LEFT)
# E1 = Entry(top, bd =5)
# E1.pack(side = RIGHT)
#
# top.mainloop()

import tkinter as tk

root = tk.Tk()
root.geometry("300x200")

def func(event):
    print("You hit return.")
root.bind('<Return>', func)

def onclick():
    print("You clicked the button")

button = tk.Button(root, text="click me", command=onclick)
button.pack()

root.mainloop()