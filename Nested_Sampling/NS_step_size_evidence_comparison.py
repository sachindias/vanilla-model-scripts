import matplotlib.pyplot as plt
import datetime

#-----------------------------------------
#MAIN EXAMPLE SCRIPT

#TIMER
start_timer = datetime.datetime.now()
print(start_timer)

rev = "0956"
spin = "afree"

step_sizes = [
            200,
            250,
            300, 
            ]

plt.figure(figsize=(7, 4))
NS_filename = "results.json"

for i in range(len(step_sizes)):
    print("%s STEPS; TIME: %s" %(step_sizes[i], datetime.datetime.now() - start_timer))
    Nested_Sampling = open("GRO_J1655_40_EPIC_rev%s_%ssteps_%s_fit/info/%s" %(rev, step_sizes[i], spin, NS_filename), "r")   
    lines = Nested_Sampling.readlines()
    Nested_Sampling.close()
    
    LnZ_val = float((lines[2].split('logz":')[1].split(","))[0])
    LnZ_err = float((lines[3].split('logzerr":')[1].split(","))[0])

    plt.errorbar(step_sizes[i], LnZ_val, yerr = LnZ_err, marker = "o", linewidth = 2, color = "k")

#PLOT SETTINGS
plt.xlabel("Number of Steps", fontsize = 14)
plt.ylabel("lnZ", fontsize = 14)

plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)

plt.title("rev%s %s - steps: %s" %(rev, spin, step_sizes), fontsize = 16)

plt.savefig("GRO_J1655_40_rev%s_%s_StepSize_Evidence_Comparison.png" %(rev, spin), bbox_inches = "tight")
plt.show()

print("TOTAL TIME: %s" %(datetime.datetime.now() - start_timer))