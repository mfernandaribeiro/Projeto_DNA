from tkinter import *

def background(seq, nl):  #nl = numero de sequencias

    janela = Tk()
    janela.title("Alinhamento Múltiplo de sequências de DNA")
    text = Text(janela)
    x = 0
    aux = ''
    while x < nl:
        text.insert(INSERT, seq[x] + '\n')
        aux = seq[x]
        for i in range(len(aux)):
            if aux[i] == 'A':
                text.tag_add("adenina" + str(x + 1), f"{x + 1}.{i}")
                text.tag_config("adenina" + str(x + 1),foreground="dark green")
            elif aux[i] == 'T':
                text.tag_add("timina" + str(x + 1), f"{x + 1}.{i}")
                text.tag_config("timina" + str(x + 1), foreground="red")
            elif aux[i] == 'C':
                text.tag_add("citosina" + str(x + 1), f"{x + 1}.{i}")
                text.tag_config("citosina" + str(x + 1), foreground="dark blue")
            elif aux[i] == 'G':
                text.tag_add("guanina" + str(x + 1), f"{x + 1}.{i}")
                text.tag_config("guanina" + str(x + 1), foreground="orange")
            elif aux[i] == '-':
                text.tag_add("gap" + str(x + 1), f"{x + 1}.{i}")
                text.tag_config("gap" + str(x + 1), foreground="purple")
            elif aux[i] == str('\n'):
                text.tag_add("n" + str(x + 1), f"{x + 1}.{i}")
                text.tag_config("n" + str(x + 1), background="black")
        x += 1
        if x == nl:
            text.insert(INSERT, seq[x])
            text.tag_add("score", f"{nl + 1}.0")
            text.pack()
    janela.mainloop()