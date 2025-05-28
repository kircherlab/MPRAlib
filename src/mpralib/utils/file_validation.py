from enum import Enum
import json
from importlib.resources import files
import jsonschema
import csv
import gzip
import ast
from mpralib.utils.io import is_compressed_file
import tqdm
import logging
import re


class ValidationSchema(Enum):
    REPORTER_SEQUENCE_DESIGN = "reporter_sequence_design"
    REPORTER_BARCODE_TO_ELEMENT_MAPPING = "reporter_barcode_to_element_mapping"
    REPORTER_EXPERIMENT_BARCODE = "reporter_experiment_barcode"
    REPORTER_EXPERIMENT = "reporter_experiment"
    REPORTER_ELEMENT = "reporter_element"
    REPORTER_VARIANT = "reporter_variant"
    REPORTER_GENOMIC_ELEMENT = "reporter_genomic_element"
    REPORTER_GENOMIC_VARIANT = "reporter_genomic_variant"


class SchemaToFileNameMap:
    def __init__(self):
        self._data = {}

    def set(self, key: ValidationSchema, file_name: str):
        if not isinstance(key, ValidationSchema):
            raise KeyError(f"Key must be a FileKey enum value. Got {key}")
        if not isinstance(file_name, str):
            raise ValueError("File path must be a string")
        self._data[key] = file_name

    def get(self, key: ValidationSchema):
        return self._data.get(key, None)

    def as_dict(self):
        return {k.value: v for k, v in self._data.items()}


schemaFilemap = SchemaToFileNameMap()
schemaFilemap.set(ValidationSchema.REPORTER_SEQUENCE_DESIGN, "reporter_sequence_design.json")
schemaFilemap.set(ValidationSchema.REPORTER_BARCODE_TO_ELEMENT_MAPPING, "reporter_barcode_to_element_mapping.json")
schemaFilemap.set(ValidationSchema.REPORTER_EXPERIMENT_BARCODE, "reporter_experiment_barcode.json")
schemaFilemap.set(ValidationSchema.REPORTER_EXPERIMENT, "reporter_experiment.json")
schemaFilemap.set(ValidationSchema.REPORTER_ELEMENT, "reporter_element.json")
schemaFilemap.set(ValidationSchema.REPORTER_VARIANT, "reporter_variant.json")
schemaFilemap.set(ValidationSchema.REPORTER_GENOMIC_ELEMENT, "reporter_genomic_element.json")
schemaFilemap.set(ValidationSchema.REPORTER_GENOMIC_VARIANT, "reporter_genomic_variant.json")


def _convert_row_value(value: str, prop_schema: dict):
    try:
        if prop_schema.get("type") == "integer":
            converted_value = int(value)
        elif prop_schema.get("type") == "number":
            converted_value = float(value)
        elif prop_schema.get("type") == "array":
            converted_value = ast.literal_eval(value)
        else:
            converted_value = value
    except ValueError:
        converted_value = value  # Let validation catch the error

    return converted_value


def validate_tsv_with_schema(tsv_file_path: str, schema_type: ValidationSchema) -> bool:
    """
    Validates a TSV file against a specified JSON schema.

    This function reads a TSV file (optionally gzipped), converts each row to a dictionary,
    and validates each row against the provided JSON schema. If any row fails validation,
    a warning is logged. If an unexpected error occurs during validation, it is logged and raised.

    Parameters:
        tsv_file_path (str): Path to the TSV file to validate. The file may be gzipped.
        schema_type (ValidationSchema): The type of schema to validate against.

    Returns:
        bool: True if all rows are valid according to the schema, False otherwise.

    Raises:
        Exception: If an unexpected error occurs during validation.

    Logs:
        - Warnings for each row that fails schema validation.
        - Errors for unexpected exceptions during validation.
        - Info if the file is valid according to the schema.
        - Warning if the file is not valid according to the schema.
    """
    LOGGER = logging.getLogger(__name__)
    LOGGER.setLevel(logging.WARNING)

    schema = _load_schema(schema_type)
    header = _get_header_for_schema(schema_type)
    open_func = gzip.open if is_compressed_file(tsv_file_path) else open
    mode = "rt" if is_compressed_file(tsv_file_path) else "r"

    correct_file = True

    with open_func(tsv_file_path, mode, encoding="utf-8") as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter="\t", fieldnames=header)  # type: ignore

        for i, row in enumerate(tqdm.tqdm(reader, desc="Validating rows", unit="row"), start=1):
            _convert_row_types(row, schema)
            try:
                jsonschema.validate(instance=row, schema=schema)
            except jsonschema.ValidationError as e:
                LOGGER.warning(f"Row {i} invalid: {e.message}")
                correct_file = False
            except Exception as e:
                LOGGER.error(f"Row {i} error: {e}")
                correct_file = False
                raise e

    if correct_file:
        LOGGER.info(f"File {tsv_file_path} is valid according to schema {schema_type.value}.")
    else:
        LOGGER.warning(f"File {tsv_file_path} is not valid according to schema {schema_type.value}.")
    return correct_file


def _load_schema(schema_type: ValidationSchema):
    schema_filename = schemaFilemap.get(schema_type)
    if schema_filename is None:
        raise ValueError(f"No schema file mapped for schema type: {schema_type}")
    schema_path = files("mpralib.schemas").joinpath(schema_filename)
    with schema_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _get_header_for_schema(schema_type: ValidationSchema) -> list | None:
    if schema_type == ValidationSchema.REPORTER_BARCODE_TO_ELEMENT_MAPPING:
        return ["barcode", "oligoName"]
    elif schema_type == ValidationSchema.REPORTER_GENOMIC_ELEMENT:
        return [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "log2FoldChange",
            "inputCount",
            "outputCount",
            "minusLog10PValue",
            "minusLog10QValue",
        ]
    elif schema_type == ValidationSchema.REPORTER_GENOMIC_VARIANT:
        return [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "log2FoldChange",
            "inputCountRef",
            "outputCountRef",
            "inputCountAlt",
            "outputCountAlt",
            "minusLog10PValue",
            "minusLog10QValue",
            "postProbEffect",
            "CI_lower_95",
            "CI_upper_95",
            "variantPos",
            "refAllele",
            "altAllele",
        ]
    return None


def _convert_row_types(row: dict, schema: dict) -> None:
    # Handle patternProperties
    for prop_pattern_string, prop_schema in schema.get("patternProperties", {}).items():
        prop_pattern = re.compile(prop_pattern_string)
        for prop in [p for p in row if prop_pattern.match(p)]:
            if row[prop] != "":
                if "anyOf" in prop_schema:
                    for anyOfProp_schema in prop_schema["anyOf"]:
                        row[prop] = _convert_row_value(row[prop], anyOfProp_schema)
                else:
                    row[prop] = _convert_row_value(row[prop], prop_schema)
    # Handle properties
    for prop, prop_schema in schema.get("properties", {}).items():
        if prop in row and row[prop] != "":
            if "anyOf" in prop_schema:
                for anyOfProp_schema in prop_schema["anyOf"]:
                    row[prop] = _convert_row_value(row[prop], anyOfProp_schema)
            else:
                row[prop] = _convert_row_value(row[prop], prop_schema)
