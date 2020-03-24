

single_count = 0
singles =[]
with open('results_single.txt','w') as out:

    with open('single_report.txt') as ins:
        for line in ins:
            if single_count==2:
                out.write(" ".join(singles)+"\n")
                single_count=0
                singles=[]
            if len(line)>3 and line[0]!='%':
                ls= line.split()
                singles.append(ls[3])
                single_count+=1