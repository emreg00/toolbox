
def main():
    file_name = "/home/emre/arastirma/data/drug/atc/br08303.keg"
    atc_to_name = get_atc_to_name(file_name)
    print atc_to_name["D07A"]
    return

def get_atc_to_name(file_name):
    atc_to_name = {}
    f = open(file_name)
    line = f.readline()
    while line[0] != "!":
        line = f.readline()
    for line in f:
	if line.startswith("!"):
	    break
	if line.startswith("#"):
	    continue
	words = line.strip().split()
        prefix = words[0]
        if len(prefix) > 1:
            code = prefix[1:]
        else:
            code = words[1]
        name = " ".join(words[2:])
        atc_to_name[code] = name
    f.close()
    return atc_to_name


if __name__ == "__main__":
    main()

