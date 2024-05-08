import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from PIL import Image, ImageTk
import webbrowser

# Define color scheme
colors = ['#FAC6D2', '#282F44', '#645070', '#9E6471', '#98A8D7']

# Function to retrieve cancer details from file
def get_cancer_details(cancer_type):
    with open('cancer_details.txt', 'r') as file:
        # Skip the header line
        next(file)

        # Parse each line in the file
        for line in file:
            # Split the line into type and description
            cancer, description = line.split('\t')

            if cancer == cancer_type:
                return description
            
    return cancer_type.lower() + "Cancer type not found"

# Function to forget all widgets in the window
def forget_all_widgets():
    # List all widgets
    all_widgets = root.winfo_children()
    # Forget all widgets
    for widget in all_widgets:
        widget.pack_forget()

# Function to display the welcome page
def show_welcome_page():
    forget_all_widgets()
    welcome_label.pack()
    website_link.pack()

    button_frame = tk.Frame(root)
    button_frame.pack(pady=10)  # Add some padding around the frame

    aggregate_button = tk.Button(button_frame, text="Aggregate", command=show_aggregate_page)
    # Dropdown menu
    cancers = ["Pancreatic", "Colon", "Lung", "Adrenal", "Blood", "Brain", "Breast", "Esophageal", "Liver", "Ovarian", "Prostate"]
    cancers.sort()
    dropdown = ttk.Combobox(button_frame, values=cancers)
    dropdown.set("Select Option")
    dropdown.bind("<<ComboboxSelected>>", lambda e: show_cancer_details(dropdown.get()))
    aggregate_button.pack(side="left", padx=10)
    dropdown.pack(side="left", padx=10)

    logo_label.pack()

# Function to display the aggregate page
def show_aggregate_page():
    forget_all_widgets()

    cancer_label.config(text="All Cancer Data", fg=colors[0], font=("Arial", 24, "bold"))
    cancer_label.pack()

    info_label.config(text="The below graphs illustrate patterns among the aggregate data collected from all cancer types")
    info_label.pack()

    # Create a frame for the images
    image_frame = tk.Frame(root)
    image_frame.pack()

    # Display two images in each row of the frame
    image_paths = ["plots/aggregate_chrom_plot.png", "plots/aggregate_mutation_type_plot.png", "plots/aggregate_transitions_transversions_plot.png", "plots/aggregate_protein_mutation_plot.png"]
    num_cols = 2
    for i in range(0, len(image_paths), num_cols):
        row_frame = tk.Frame(image_frame)
        row_frame.pack()
        for j in range(num_cols):
            if i + j < len(image_paths):
                img = Image.open(image_paths[i + j])
                img = img.resize((350, 300))
                photo = ImageTk.PhotoImage(img)
                image_label = tk.Label(row_frame, image=photo)
                image_label.image = photo  # Prevent garbage collection
                image_label.pack(side="left", padx=10, pady=10)

    home_button.pack()

# Function to display the image
def show_image(prot_or_chrom_or_tran, cancer_type):
    info_label.pack_forget()
    image_path = "plots/" + cancer_type.lower() + prot_or_chrom_or_tran + "_plot.png"
    
    # Load image
    image = Image.open(image_path)
    photo = ImageTk.PhotoImage(image)
    
    # Create a label to display the image
    image_label.config(image=photo)
    image_label.image = photo
    image_label.pack()

# Function to display cancer details
def show_cancer_details(cancer_type):
    forget_all_widgets()
    cancer_label.config(text=f"{cancer_type} Cancer", fg=colors[0], font=("Arial", 24, "bold"))
    cancer_label.pack()
    info_label.config(text=f"{get_cancer_details(cancer_type.lower())}", font=("Arial", 18))
    info_label.pack()

    # Buttons
    button_frame = tk.Frame(root)
    button_frame.pack(pady=10)  # Add some padding around the frame
    protein_button = tk.Button(button_frame, text="Protein Mutations", command=lambda: show_image("_protein_mutation", cancer_type))
    chromosome_button = tk.Button(button_frame, text="Chromosomal Mutations", command=lambda: show_image("_chrom", cancer_type))
    tran_button = tk.Button(button_frame, text="Transitions vs Transversions", command=lambda: show_image("_transitions_transversions", cancer_type))
    mutation_type_button = tk.Button(button_frame, text="Mutation Types", command=lambda: show_image("_mutation_type", cancer_type))
    protein_button.pack(side="left", padx=10)
    chromosome_button.pack(side="left", padx=10)
    tran_button.pack(side="left", padx=10)
    mutation_type_button.pack(side="left", padx=10)
    home_button.pack()

# Create the main window
root = tk.Tk()
root.title("Welcome Page")

# Welcome label and website link
welcome_label = tk.Label(root, text="Welcome!", fg=colors[0], font=("Arial", 50, "bold"))
website_link = tk.Label(root, text="Data Collected from the GDC Data Portal", fg=colors[4], cursor="hand2")
website_link.bind("<Button-1>", lambda e: webbrowser.open("https://portal.gdc.cancer.gov/"))

# Buttons
home_button = tk.Button(root, text="Home", command=show_welcome_page)

# Cancer type label
cancer_label = tk.Label(root)

# Information label
info_label = tk.Label(root, wraplength=600)

# Image for specific image
image_label = tk.Label(root)

logo = Image.open("data/gdc_logo.png")
photo = ImageTk.PhotoImage(logo)

# Create a label to display the image
logo_label = tk.Label(root, image=photo)

# Show the welcome page
show_welcome_page()

# Setting window size
root.geometry("800x800")

# Run the application
root.mainloop()
