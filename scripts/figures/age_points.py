"""Read and verify the patient-level age coordinates for Figure 1B."""
import re

from openpyxl import load_workbook
from scipy.stats import spearmanr

AGE_FIELDS = {
    "context": "Clinical context age (years)",
    "first": "First event age (years)",
    "hematologic": "Hematologic event age (years)",
    "mixed": "Mixed reported age (years)",
}


def read_clinical(workbook_path):
    workbook = load_workbook(workbook_path, data_only=True, read_only=True)
    cohort = workbook["S1 Cohort"]
    values = list(cohort.values)
    headers = values[4]
    patients = [dict(zip(headers, row)) for row in values[5:183]]

    # Broader reported-age records ineligible for the clinical-context plot
    # still have a primary germline residue in the S2A allele links.
    primary_position = {}
    for row in list(workbook["S2 Variants"].values)[5:190]:
        patient, _, protein = row[:3]
        if patient and protein:
            primary_position.setdefault(patient, int(re.search(r"\d+", protein).group()))
    workbook.close()
    groups = {}
    for endpoint, field in AGE_FIELDS.items():
        eligible = "Mixed age-position eligible" if endpoint == "mixed" else "Age-position eligible"
        coords = []
        for row in patients:
            age = row[field]
            if row[eligible] != "Yes" or age is None:
                continue
            patient = row["Patient ID"]
            position = row["Amino-acid position"]
            if position is None:
                position = primary_position[patient]
            coords.append((int(position), float(age), patient))
        groups[endpoint] = coords
    return groups


def rank_summary(points):
    positions, ages, _ = zip(*points)
    result = spearmanr(positions, ages)
    return {"n": len(points), "rho": float(result.statistic), "p": float(result.pvalue)}


