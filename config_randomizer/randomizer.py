from string import Template 
import random 


# TODO discuss/add randomized parameters/ranges


replacements = {
    'mt_quantity': random.randint(1, 100),
    'mt_length': random.randint(1, 100),
    'linker_quantity': random.randint(1,100),
    'blob_quantity': random.randint(1,50),


}

with open('configtemplate', 'r') as file:
    src = Template(file.read())
    print(src)
    result = src.substitute(replacements)
    print(result)
    config = open('config.cym', 'w')
    config.write(result)
    config.close()
file.close()

