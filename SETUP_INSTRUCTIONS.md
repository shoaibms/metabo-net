# **📌 SETUP_INSTRUCTIONS.md**
> **Complete Setup Guide for the Data Cleaning Pipeline**

---

## **🔹 Option 1: Install using `pip` (Recommended for Python Users)**
If you are using **Python only**, you can install all dependencies with:
```bash
pip install -r requirements.txt
```

### **✔ Verify Installation**
After installation, you can verify if all required packages are installed:
```bash
python -c "import pandas, numpy, matplotlib, seaborn, sklearn, scipy; print('All packages installed successfully!')"
```

---
## **🔹 Option 2: Install using Conda (Recommended for Python & R Users)**
If you want **both Python and R dependencies**, use Conda for a fully controlled environment.

### **✔ Step 1: Create a Conda Environment**
```bash
conda env create -f environment.yml
```

### **✔ Step 2: Activate the Environment**
```bash
conda activate data-cleaning-pipeline
```

### **✔ Step 3: Verify Installation**
Run:
```bash
conda list
```
This should list all installed dependencies.

---
## **🔹 Option 3: Install R Dependencies Manually**
If you're using **R for imputation**, install the required packages manually by opening **R** and running:

```r
install.packages(c("dplyr", "missForest", "doParallel", "mice"))
```

---

## **🔹 Option 4: Using Docker (Optional)**
If you prefer **Docker** for full reproducibility, follow these steps:

### **✔ Step 1: Build the Docker Image**
```bash
docker build -t data-cleaning-pipeline .
```

### **✔ Step 2: Run the Container**
```bash
docker run --rm -v $(pwd):/app data-cleaning-pipeline
```

This will run the pipeline inside a container with all dependencies installed.

---

## **🔹 Frequently Asked Questions**
### **1️⃣ Which setup should I use?**
- **If you only need Python**, use **pip** (`requirements.txt`).
- **If you need both Python & R**, use **Conda** (`environment.yml`).
- **If you need full isolation & reproducibility**, use **Docker**.

### **2️⃣ What if a package is missing?**
Run:
```bash
pip install <missing_package>
```
Or, if using Conda:
```bash
conda install <missing_package>
```

### **3️⃣ How do I remove the Conda environment?**
If you want to delete the Conda environment, run:
```bash
conda env remove -n data-cleaning-pipeline
```

---

## **🎯 Final Notes**
✅ Choose the best setup based on your needs.  
✅ Follow the installation steps carefully.  
✅ If you face issues, check **Python/R versions** and ensure **internet access**.  

🚀 **Happy Data Cleaning!** 🎉
