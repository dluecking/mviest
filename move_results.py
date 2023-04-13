import glob
import os.path
import shutil

samples = glob.glob("results/*")
samples = [x.replace("results/", "") for x in samples]

results_done = []
results_not_done = []

for sample in samples:
    png_file = "results/" + sample + "/" + sample + "_mviest_plot.png"

    if os.path.isfile(png_file):
        results_done.append(sample)
    else:
        results_not_done.append(sample)

print("Not finished:")
print(results_not_done)

print("Done:")
print(results_done)

move = input("Should we move the finished results to ~/bioinf/mviest_results/tara ? <yes/no> \n")

if move == "yes":
    for sample in results_done:
        src = "results/" + sample
        dest = "/bioinf/home/dlueckin/mviest_results/tara"
        shutil.move(src, dest)
