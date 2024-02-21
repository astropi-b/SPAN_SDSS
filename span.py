import tkinter as tk
from tkinter import filedialog, Toplevel, simpledialog, messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from astropy.io import fits
import numpy as np
from scipy.constants import k, m_e
from scipy.signal import find_peaks
import tempfile
import os
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import textwrap



# Constants for physics calculations
b = 2.897e-3  # Wien's displacement constant in m*K

def moving_average(data, window_size):
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(data, window, 'same')

class SpectrumAnalyzer(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Spectrum Analyzer: Group Quazars")
        self.geometry("1200x1200")
        self.create_widgets()
        self.pdf_path = simpledialog.askstring("PDF Path", "Enter the path for the PDF report:", parent=self)
        # If the user does not enter a path, use a default value
        if not self.pdf_path:
            self.pdf_path = "spectrum_analysis_report.pdf"
        self.plots_paths = []  # To store paths of plot images
        self.report_items = []
        # Initialize the analysis_texts list to hold text items for the report
        self.analysis_texts = []
        # Example of how you might populate this attribute
        self.populate_analysis_texts()


    def populate_analysis_texts(self):
        # Just an example - you would replace this with your actual method of populating texts
        self.analysis_texts.append("Results.")

    def show_analysis(self):
        # Example method that uses the analysis_texts attribute
        for text in self.analysis_texts:
            print(text)  # or display it in your GUI


    def create_widgets(self):
        self.upload_button = tk.Button(self, text="Upload FITS File", command=self.upload_and_plot_spectrum)
        self.upload_button.pack(side=tk.TOP, pady=20)
        self.generate_pdf_button = tk.Button(self, text="Generate PDF Report", command=self.generate_pdf_report)
        self.generate_pdf_button.pack(side=tk.TOP, pady=10)

    def upload_and_plot_spectrum(self):
        file_path = filedialog.askopenfilename(title="Select the FITS file of the star spectrum")
        if not file_path:
            return

        with fits.open(file_path) as hdul:
            data = hdul[1].data
            flux = data['flux']
            loglam = data['loglam']
            wavelength = 10 ** loglam

        self.plot_spectrum(wavelength, flux, 'Spectral Data', 'Original Spectrum')
        self.ask_window_size_and_smooth(wavelength, flux)

    def plot_spectrum(self, wavelength, flux, label, title, color='blue'):
        top = Toplevel(self)
        top.title(title)
        fig = Figure(figsize=(20, 5), dpi=100)
        plot = fig.add_subplot(1, 1, 1)
        plot.plot(wavelength, flux, label=label, color=color)
        plot.set_xlabel('Wavelength (Angstrom)')
        plot.set_ylabel('Flux')
        plot.set_title(title)
        plot.grid(True)
        plot.legend()

        canvas = FigureCanvasTkAgg(fig, master=top)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Save the figure for the PDF report
        plot_img_path = tempfile.mktemp(suffix=".png")
        fig.savefig(plot_img_path)
        self.plots_paths.append(plot_img_path)
    
    def ask_window_size_and_smooth(self, wavelength, flux):
        window_size = simpledialog.askinteger("Input", "Enter the window size for smoothing:", parent=self)
        if window_size is not None:
            smoothed_flux = moving_average(flux, window_size)
            
            self.plot_avgspectrum(wavelength, smoothed_flux)
            self.plot_spectrum_with_absorption_lines(wavelength, smoothed_flux)
            self.plot_spectrum_with_emission_lines(wavelength, smoothed_flux)
            self.calculate_temperature_and_plot_mb_distribution(wavelength, smoothed_flux)
            self.plot_spectrum_with_absspec_lines(wavelength, smoothed_flux)
            self.plot_spectrum_with_emispec_lines(wavelength, smoothed_flux)
            self.saha_equation_h(wavelength,smoothed_flux)
            self.saha_equation_ca(wavelength,smoothed_flux)
            self.boltzmann_equation_h(wavelength,smoothed_flux)
            self.boltzmann_equation_ca(wavelength,smoothed_flux)
            

    def plot_avgspectrum(self, wavelength, smoothed_flux):
        inverted_smoothed_flux = -smoothed_flux
        peaks, _ = find_peaks(inverted_smoothed_flux, prominence=1.55)  # Adjust prominence as needed

        # Plot the spectrum and mark the absorption lines
        top = Toplevel(self)
        top.title("Star Spectrum Averaged")
        fig = Figure(figsize=(10, 4), dpi=100)
        plot = fig.add_subplot(1, 1, 1)
        plot.plot(wavelength, smoothed_flux, label='Smoothed Spectral Data') 
        plot.set_xlabel('Wavelength (Angstrom)')
        plot.set_ylabel('Flux')
        plot.set_title('Star Avg Spectrum')
        plot.legend()
        plot.grid(True)

        canvas = FigureCanvasTkAgg(fig, master=top)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Save the figure for the PDF report
        plot_img_path = tempfile.mktemp(suffix=".png")
        fig.savefig(plot_img_path)
        self.plots_paths.append(plot_img_path)
    
    def plot_spectrum_with_absorption_lines(self, wavelength, smoothed_flux):
        # Find absorption lines by finding peaks in the inverted flux
        inverted_smoothed_flux = -smoothed_flux
        peaks, _ = find_peaks(inverted_smoothed_flux, prominence=1.55)  # Adjust prominence as needed

        # Plot the spectrum and mark the absorption lines
        top = Toplevel(self)
        top.title("Star Spectrum with Absorption Lines")
        fig = Figure(figsize=(10, 4), dpi=100)
        plot = fig.add_subplot(1, 1, 1)
        plot.plot(wavelength, smoothed_flux, label='Smoothed Spectral Data')
        plot.scatter(wavelength[peaks], smoothed_flux[peaks], color='red', s=40, label='Absorption Lines', zorder=5)
        plot.set_xlabel('Wavelength (Angstrom)')
        plot.set_ylabel('Flux')
        plot.set_title('Star Spectrum with Absorption Lines')
        plot.legend()
        plot.grid(True)

        canvas = FigureCanvasTkAgg(fig, master=top)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Save the figure for the PDF report
        plot_img_path = tempfile.mktemp(suffix=".png")
        fig.savefig(plot_img_path)
        self.plots_paths.append(plot_img_path)

        # Display wavelengths and fluxes of absorption lines
        # Limit to the first 10 peaks for display
        absorption_lines_info = "\n".join([f"Wavelength: {wavelength[p]:.2f} Å, Flux: {smoothed_flux[p]:.2f}" for p in peaks])
        max_display = 10
        limited_absorption_lines_info = "\n".join([f"Wavelength: {wavelength[p]:.2f} Å, Flux: {smoothed_flux[p]:.2f}" for p in peaks[:max_display]])

        # Check if there are more lines than displayed
        if len(peaks) > max_display:
            limited_absorption_lines_info += "\n...and more lines not shown here."

        messagebox.showinfo("Absorption Lines", f"Absorption Lines Detected:\n{limited_absorption_lines_info}")

        
        #messagebox.showinfo("Absorption Lines", f"Absorption Lines Detected:\n{absorption_lines_info}")
        # Example of adding formatted text for wavelengths and fluxes
        # Assuming 'peaks' is an array of indices where the peaks occur in the wavelength and smoothed_flux arrays
        for peak_index in peaks:
            entry = f"Absorption line Wavelength: {wavelength[peak_index]:.2f} Å, Absorption line Flux: {smoothed_flux[peak_index]:.2f}"
            self.analysis_texts.append(entry)

    def plot_spectrum_with_absspec_lines(self, wavelength, smoothed_flux):
        # Find absorption lines by finding peaks in the inverted flux
        inverted_smoothed_flux = -smoothed_flux
        peaks, _ = find_peaks(inverted_smoothed_flux, prominence=1.55)# Adjust prominence as needed
        known_lines = {
        "Mg I": 5175,
        "Hα": 6562.8,
        "Hβ": 4861.3,
        "H gamma": 4340.5,
        "Ca II (K)": 3933.7,
        "Ca II (H)": 3968.5,
        "He I (5876)": 5875.6,
        "He I (6678)": 6678.2,
        "K I (7665)": 7664.9,
        "K I (7700)": 7699.0,
        "Na I (D1)": 5895.6,  # Sodium D1 line
        "Na I (D2)": 5889.95,  # Sodium D2 line
        "O I": 7773,  # Oxygen triplet around 7774 Å
         "Si II": 6347.1,  # Silicon, a common line in hotter stars
         "Fe II": 5169,  # Iron, indicative of various physical processes
        "N II": 6584,  # Nitrogen, seen in emission in some nebulae or active regions
        "S II": 6716,  # Sulfur, often seen in emission in nebulae
        "Fe I": 5270,  # Neutral iron, common in cooler stars
        "Ti II": 4501  # Titanium, seen in stars with rich metal lines
        }
        # Define a region width in Angstroms for highlighting
        region_width = 2.5
        tolerance = 7.0 
        # Plot the spectrum and mark the absorption lines
        top = Toplevel(self)
        top.title("Star Spectrum with Spectral Lines of elements")
        fig = Figure(figsize=(10, 4), dpi=100)
        plot = fig.add_subplot(1, 1, 1)
        plot.plot(wavelength, smoothed_flux, label='Smoothed Spectral Data')
        plot.scatter(wavelength[peaks], smoothed_flux[peaks], color='red', s=40, label='Absorption Lines', zorder=5)
        # Highlight regions and label the lines
        identified_lines = []
        for line_name, line_wavelength in known_lines.items():
    # Find if there's a peak near the line within the defined tolerance
            nearby_peak = any(abs(line_wavelength - wavelength[peak]) <= tolerance for peak in peaks)
            if nearby_peak:
                identified_lines.append(line_name)
                region_start = line_wavelength - region_width / 2
                region_end = line_wavelength + region_width / 2

                # Draw a shaded region around the line
                plot.axvspan(region_start, region_end, color='purple', alpha=0.5)
        
        # Label the line
                plot.text(line_wavelength, max(smoothed_flux), line_name, horizontalalignment='center', verticalalignment='bottom')
        
        
        plot.set_xlabel('Wavelength (Angstrom)')
        plot.set_ylabel('Flux')
        plot.set_title('Star Spectrum with Highlighted Absorption Lines')
        plot.legend()
        plot.grid(True)

        canvas = FigureCanvasTkAgg(fig, master=top)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        if identified_lines:
            identified_lines_str = "\n".join(identified_lines)
            messagebox.showinfo("Identified Absorption Lines", f"Identified absorption lines in the spectrum:\n{identified_lines_str}")
            self.analysis_texts.append(f"Identified Absorption Lines:\n{identified_lines_str}")

        else:
            messagebox.showinfo("Identified Absorption Lines", "No absorption lines identified in the spectrum.")
            self.analysis_texts.append("No absorption lines identified in the spectrum.")

        # Save the figure for the PDF report
        plot_img_path = tempfile.mktemp(suffix=".png")
        fig.savefig(plot_img_path)
        self.plots_paths.append(plot_img_path)

    def plot_spectrum_with_emispec_lines(self, wavelength, smoothed_flux):
        peaks1, _ = find_peaks(smoothed_flux, prominence=1.55)# Adjust prominence as needed
        known_lines = {
        "Mg I": 5175,
        "Hα": 6562.8,
        "Hβ": 4861.3,
        "H gamma": 4340.5,
        "Ca II (K)": 3933.7,
        "Ca II (H)": 3968.5,
        "He I (5876)": 5875.6,
        "He I (6678)": 6678.2,
        "K I (7665)": 7664.9,
        "K I (7700)": 7699.0,
        "Na I (D1)": 5895.6,  # Sodium D1 line
        "Na I (D2)": 5889.95,  # Sodium D2 line
        "O I": 7773,  # Oxygen triplet around 7774 Å
         "Si II": 6347.1,  # Silicon, a common line in hotter stars
         "Fe II": 5169,  # Iron, indicative of various physical processes
        "N II": 6584,  # Nitrogen, seen in emission in some nebulae or active regions
        "S II": 6716,  # Sulfur, often seen in emission in nebulae
        "Fe I": 5270,  # Neutral iron, common in cooler stars
        "Ti II": 4501  # Titanium, seen in stars with rich metal lines
        }
        # Define a region width in Angstroms for highlighting
        region_width = 2.5
        tolerance = 7.0 
        # Plot the spectrum and mark the absorption lines
        top = Toplevel(self)
        top.title("Star Spectrum with Emsion Spectral Lines of elements")
        fig = Figure(figsize=(10, 4), dpi=100)
        plot = fig.add_subplot(1, 1, 1)
        plot.plot(wavelength, smoothed_flux, label='Smoothed Spectral Data')
        plot.scatter(wavelength[peaks1], smoothed_flux[peaks1], color='red', s=40, label='Emission Lines', zorder=5)
        # Highlight regions and label the lines
        identified_lines = []
        for line_name, line_wavelength in known_lines.items():
    # Find if there's a peak near the line within the defined tolerance
            nearby_peak = any(abs(line_wavelength - wavelength[peak]) <= tolerance for peak in peaks1)
            if nearby_peak:
                identified_lines.append(line_name)
                region_start = line_wavelength - region_width / 2
                region_end = line_wavelength + region_width / 2

                # Draw a shaded region around the line
                plot.axvspan(region_start, region_end, color='purple', alpha=0.5)
        
        # Label the line
                plot.text(line_wavelength, max(smoothed_flux), line_name, horizontalalignment='center', verticalalignment='bottom')
        
        
        plot.set_xlabel('Wavelength (Angstrom)')
        plot.set_ylabel('Flux')
        plot.set_title('Star Spectrum with Highlighted Emission Lines')
        plot.legend()
        plot.grid(True)

        canvas = FigureCanvasTkAgg(fig, master=top)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        if identified_lines:
            identified_lines_str = "\n".join(identified_lines)
            messagebox.showinfo("Identified Emission Lines", f"Identified emission lines in the spectrum:\n{identified_lines_str}")
            self.analysis_texts.append(f"Identified Emission Lines:\n{identified_lines_str}")

        else:
            messagebox.showinfo("Identified Emission Lines", "No absorption lines identified in the spectrum.")
            self.analysis_texts.append("No emission lines identified in the spectrum.")

        # Save the figure for the PDF report
        plot_img_path = tempfile.mktemp(suffix=".png")
        fig.savefig(plot_img_path)
        self.plots_paths.append(plot_img_path)


    def plot_spectrum_with_emission_lines(self, wavelength, smoothed_flux):
        # Find absorption lines by finding peaks in the inverted flux
        peaks1, _ = find_peaks(smoothed_flux, prominence=1.55)

        # Plot the spectrum and mark the absorption lines
        top = Toplevel(self)
        top.title("Star Spectrum with Emission Lines")
        fig = Figure(figsize=(10, 4), dpi=100)
        plot = fig.add_subplot(1, 1, 1)
        plot.plot(wavelength, smoothed_flux, label='Smoothed Spectral Data')
        plot.scatter(wavelength[peaks1], smoothed_flux[peaks1], color='red', s=40, label='Emission Lines', zorder=5)
        plot.set_xlabel('Wavelength (Angstrom)')
        plot.set_ylabel('Flux')
        plot.set_title('Star Spectrum with Emission Lines')
        plot.legend()
        plot.grid(True)

        canvas = FigureCanvasTkAgg(fig, master=top)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Save the figure for the PDF report
        plot_img_path = tempfile.mktemp(suffix=".png")
        fig.savefig(plot_img_path)
        self.plots_paths.append(plot_img_path)

        # Display wavelengths and fluxes of absorption lines
        emission_lines_info = "\n".join([f"Wavelength: {wavelength[p]:.2f} Å, Flux: {smoothed_flux[p]:.2f}" for p in peaks1])
        #messagebox.showinfo("Emission Lines", f"Emission Lines Detected:\n{emission_lines_info}")
        max_display = 10
        limited_emission_lines_info = "\n".join([f"Wavelength: {wavelength[p]:.2f} Å, Flux: {smoothed_flux[p]:.2f}" for p in peaks1[:max_display]])

        # Check if there are more lines than displayed
        if len(peaks1) > max_display:
            limited_emission_lines_info += "\n...and more lines not shown here."

        messagebox.showinfo("Emission Lines", f"Emission Lines Detected:\n{limited_emission_lines_info}")
        # Example of adding formatted text for wavelengths and fluxes
        # Assuming 'peaks' is an array of indices where the peaks occur in the wavelength and smoothed_flux arrays
        for peak_index in peaks1:
            entry = f"Emission Peak Wavelength: {wavelength[peak_index]:.2f} Å, Emission Peak Flux: {smoothed_flux[peak_index]:.2f}"
            self.analysis_texts.append(entry)



    def calculate_temperature_and_plot_mb_distribution(self, wavelength, smoothed_flux):
        # Convert wavelength to meters for Wien's law
        wavelength_meters = wavelength * 1e-10

        # Find the wavelength of the maximum flux
        max_flux_index = np.argmax(smoothed_flux)
        lambda_max = wavelength_meters[max_flux_index]

        # Calculate the temperature using Wien's displacement law
        temperature = b / lambda_max

        # Display the estimated temperature using messagebox
        messagebox.showinfo("Estimated Temperature", f"Estimated temperature of the star: {temperature:.2f} K")
        self.analysis_texts.append(f"Estimated Temperature of Star is : {temperature:.2f}")
        # Plot the Maxwell-Boltzmann Distribution for the estimated temperature
        self.plot_maxwell_boltzmann_distribution(temperature)
 
    def plot_maxwell_boltzmann_distribution(self, T):
        v = np.linspace(0, 2e6, 1000)  # Speed range
        distribution = 4 * np.pi * (m_e / (2 * np.pi * k * T))**(3/2) * v**2 * np.exp(-m_e * v**2 / (2 * k * T))
        Vmp = np.sqrt(2 * k * T / m_e)
        Vrms = np.sqrt(3 * k * T / m_e)
        top = Toplevel(self)
        top.title("Maxwell-Boltzmann Distribution")
        fig = Figure(figsize=(8, 3), dpi=100)
        plot = fig.add_subplot(1, 1, 1)
        plot.plot(v, distribution, label='Distribution')
        plot.axvline(x=Vmp, color='r', linestyle='--', label=f'Vmp = {Vmp:.2e} m/s')
        plot.axvline(x=Vrms, color='g', linestyle='--', label=f'Vrms = {Vrms:.2e} m/s')
        plot.set_xlabel('Speed (m/s)')
        plot.set_ylabel('Probability density')
        plot.grid(True)
        plot.set_title('Maxwell-Boltzmann Distribution')

        canvas = FigureCanvasTkAgg(fig, master=top)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Save the figure for the PDF report
        plot_img_path = tempfile.mktemp(suffix=".png")
        fig.savefig(plot_img_path)
        self.plots_paths.append(plot_img_path)
        messagebox.showinfo("Most Probable Speed ", f"Estimated Vmp: {Vmp:.2e} m/s")
        messagebox.showinfo("Root-Mean-Square Speed ", f"Estimated Vrms: {Vrms:.2e} m/s")
        self.analysis_texts.append(f"Most Probable Speed is : {Vmp:.2e} m/s")
        self.analysis_texts.append(f"Root-Mean-Square Speed is : {Vrms:.2e} m/s")

    def saha_equation_h(self,wavelength,smoothed_flux):
        wavelength_meters = wavelength * 1e-10

        # Find the wavelength of the maximum flux
        max_flux_index = np.argmax(smoothed_flux)
        lambda_max = wavelength_meters[max_flux_index]

        # Calculate the temperature using Wien's displacement law
        temperature = b / lambda_max
        T=temperature
        electron_pressure_str = simpledialog.askstring("Input", "Enter the electron pressure (e.g., 1e5, 1e-14):", parent=self)
        try:
            # Convert the string input to a floating-point number
            electron_pressure = float(electron_pressure_str)
            print("Electron Pressure:", electron_pressure)
        except ValueError:
            # Handle the case where the input is not a valid float
            print("Invalid input. Please enter a number in scientific notation (e.g., 1e5, 1e-14).")

        k = 1.3806452e-23  # Boltzmann constant in eV/K
        h = 6.62607015e-34  # Planck constant
        m_e = 9.10938356e-31  # Electron mass
        ionization_energy = 13.6
    
        part1 = (2 * np.pi * m_e * k * T) / (h**2)
        part2 = part1**(3/2)
        part3 = (2 *k*T* 1) / 2
        part4 = np.exp((-ionization_energy)* (1.60218e-19)/ (k * T))
        
        saha_factor = (part3 * part4 * part2) / electron_pressure
    
        
        messagebox.showinfo("Ratio of ionized to neutral hydrogen ", f"Ratio of ionized to neutral hydrogen : {saha_factor}")
        self.analysis_texts.append(f"Ratio of ionized to neutral hydrogen : {saha_factor}")

    def saha_equation_ca(self,wavelength,smoothed_flux):
        wavelength_meters = wavelength * 1e-10

        # Find the wavelength of the maximum flux
        max_flux_index = np.argmax(smoothed_flux)
        lambda_max = wavelength_meters[max_flux_index]

        # Calculate the temperature using Wien's displacement law
        temperature = b / lambda_max
        T=temperature
        electron_pressure_str = simpledialog.askstring("Input", "Enter the electron pressure (e.g., 1e5, 1e-14):", parent=self)
        try:
            # Convert the string input to a floating-point number
            electron_pressure = float(electron_pressure_str)
            print("Electron Pressure:", electron_pressure)
        except ValueError:
            # Handle the case where the input is not a valid float
            print("Invalid input. Please enter a number in scientific notation (e.g., 1e5, 1e-14).")

        k = 1.3806452e-23  # Boltzmann constant in eV/K
        h = 6.62607015e-34  # Planck constant
        m_e = 9.10938356e-31  # Electron mass
        ionization_energy = 6.11 
    
        part1 = (2 * np.pi * m_e * k * T) / (h**2)
        part2 = part1**(3/2)
        part3 = (2 *k*T* 2.30) / 1.32
        part4 = np.exp((-ionization_energy)* (1.60218e-19)/ (k * T))
    
        saha_factor = (part3 * part4 * part2) / electron_pressure
        
        messagebox.showinfo("Ratio of ionized to Calcium ", f"Ca Ionization Ratio: {saha_factor}")
        self.analysis_texts.append(f"Ca Ionization Ratio : {saha_factor}")


    def boltzmann_equation_h(self,wavelength,smoothed_flux):
        wavelength_meters = wavelength * 1e-10

        # Find the wavelength of the maximum flux
        max_flux_index = np.argmax(smoothed_flux)
        lambda_max = wavelength_meters[max_flux_index]

        # Calculate the temperature using Wien's displacement law
        temperature = b / lambda_max
        T=temperature
        k = 1.3806452e-23  # Boltzmann constant
        h = 6.62607015e-34  # Planck constant
        m_e = 9.10938356e-31  # Electron mass
        energy_difference=10.2
        part1 = 8 / 2
        part2 = np.exp(-(energy_difference)* (1.60218e-19) / (k * T))
        boltzman_factor = (part1*part2)
        messagebox.showinfo("H Excitation Ratio: ", f"H Excitation Ratio: {boltzman_factor}")
        self.analysis_texts.append(f"H Excitation Ratio: {boltzman_factor}")

    def boltzmann_equation_ca(self,wavelength,smoothed_flux):
        wavelength_meters = wavelength * 1e-10

        # Find the wavelength of the maximum flux
        max_flux_index = np.argmax(smoothed_flux)
        lambda_max = wavelength_meters[max_flux_index]

        # Calculate the temperature using Wien's displacement law
        temperature = b / lambda_max
        T=temperature
        k = 1.3806452e-23  # Boltzmann constant
        h = 6.62607015e-34  # Planck constant
        m_e = 9.10938356e-31  # Electron mass
        energy_difference=0
        part1 = 2 / 2.32
        part2 = np.exp((energy_difference)* (1.60218e-19) / (k * T))
        boltzman_factor = (part1*part2)
        messagebox.showinfo("CaI Excitation Ratio: ", f"Ca Excitation Ratio: {boltzman_factor}")
        self.analysis_texts.append(f"CaI Excitation Ratio: {boltzman_factor}")
        
        
        
        

    def generate_pdf_report(self):
        c = canvas.Canvas(self.pdf_path, pagesize=letter)
        width, height = letter

        for plot_path in self.plots_paths:
            c.drawImage(plot_path, 50, height - 250, width=400, height=200, preserveAspectRatio=True)
            height -= 250  # Adjust spacing as needed
            if height < 250:
                c.showPage()
                height = letter[1]
        # Add analysis texts to the PDF with improved formatting
        text_margin = 15  # Adjust based on your preference for spacing between lines
        for text_item in self.analysis_texts:
            text_str = str(text_item)
        # Split text into multiple lines if it's too long
            wrapped_text = textwrap.wrap(text_str, width=80)  # Adjust width as needed
            for line in wrapped_text:
                c.drawString(50, height - 20, line)
                height -= text_margin  # Move to the next line
                if height < 100:  # Check if we need a new page
                    c.showPage()
                    height = letter[1]

        c.save()
        messagebox.showinfo("PDF Report", "The PDF report has been generated.")

# Update your methods that generate analysis results to append text results to self.analysis_texts for PDF reporting

if __name__ == "__main__":
    app = SpectrumAnalyzer()
    app.mainloop()
