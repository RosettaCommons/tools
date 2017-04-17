data = []
for line in open('ST.job_list'):
    if '150000' in line:
        data.append((line, 0))
    else:
        data.append((line, 1))

data = sorted(data, key=lambda x: x[1])
out = open('ST.job_list_sorted', 'w')
for a in data:
    out.write(a[0])
out.close()
