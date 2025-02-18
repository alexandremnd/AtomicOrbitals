import numpy as np
import os
from atom_list import get_element_charge, get_element_name, atom_list

def get_file_content(file_path) -> list[str]:
    with open(file_path, 'r') as file:
        return file.readlines()

def get_beginning(file_content: list[str]) -> int:
    for i, line in enumerate(file_content):
        splitted_line = line.split()
        if (len(splitted_line) > 0 and get_element_charge(splitted_line[0]) != -1):
            return i
    return -1

def parse_atom(atom_lines: list[str]) -> dict:
    atom = {
        'name': get_element_name(atom_lines[0].split()[0]),
        'charge': get_element_charge(atom_lines[0].split()[0]),
        'basis': []
    }

    current_orbital_type = ""
    current_orbital_coefficients = []
    current_orbital_coefficients_bis = []
    current_orbital_exponents = []

    for line_no in range(1, len(atom_lines)):
        splitted = atom_lines[line_no].split()

        if len(splitted) == 0:
            continue

        if splitted[0] in ["S", "P", "D", "F", "SP"]:
            if current_orbital_type != "" and current_orbital_type != "SP":
                atom['basis'].append((current_orbital_type, current_orbital_coefficients, current_orbital_exponents))

            if current_orbital_type != "" and current_orbital_type == "SP":
                atom['basis'].append(("S", current_orbital_coefficients, current_orbital_exponents))
                atom['basis'].append(("P", current_orbital_coefficients_bis, current_orbital_exponents))

            current_orbital_type = splitted[0]
            current_orbital_coefficients = []
            current_orbital_coefficients_bis = []
            current_orbital_exponents = []
        else:
            current_orbital_exponents.append(float(splitted[0].replace("D", "E")))
            current_orbital_coefficients.append(float(splitted[1].replace("D", "E")))

            if current_orbital_type == "SP":
                current_orbital_coefficients_bis.append(float(splitted[1].replace("D", "E")))

    # Add the last orbital
    if current_orbital_type != "" and current_orbital_type != "SP":
        atom['basis'].append((current_orbital_type, current_orbital_coefficients, current_orbital_exponents))

    if current_orbital_type != "" and current_orbital_type == "SP":
        atom['basis'].append(("S", current_orbital_coefficients, current_orbital_exponents))
        atom['basis'].append(("P", current_orbital_coefficients_bis, current_orbital_exponents))

    print(atom)
    return atom


def parse_file(file_path):
    # WARNING: This function is expecting Gaussian file type.
    file_content: list[str] = get_file_content(file_path)
    first_line = get_beginning(file_content)
    file_str: str = ''.join(file_content[first_line:-1])

    atom_declaration = file_str.split('****\n')
    atoms = []

    print(f"Parsing {len(atom_declaration)} atoms from {file_path}")

    for declaration in atom_declaration:
        atoms.append(parse_atom(declaration.split("\n")))

    return atoms

def write_basis_file(atom: dict, basis_name: str) -> str:
    basis_file_name = f"{atom['name'].lower()}.basis"
    basis_file = ""

    for orbital in atom['basis']:
        basis_file += f"{orbital[0]} {len(orbital[1])} {orbital[1]} {orbital[2]} \n".replace("[", "").replace("]", "").replace(",", "")

    if not os.path.exists(f"./../data/basis_set/{basis_name}"):
        os.makedirs(f"./../data/basis_set/{basis_name}")

    with open(f"./../data/basis_set/{basis_name}/{basis_file_name}", "w+") as file:
        file.write(basis_file)


basis_name = "sto6g"
atoms = parse_file(f"./{basis_name}.txt")

for atom in atoms:
    write_basis_file(atom, basis_name)