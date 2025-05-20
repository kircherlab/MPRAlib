from enum import Enum
import json
from importlib.resources import files
import jsonschema
import csv
import gzip
import ast
from mpralib.utils.io import is_compressed_file
import tqdm


class ValidationSchema(Enum):
    REPORTER_SEQUENCE_DESIGN = "reporter_sequence_design"


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


def _convert_row_value(value: str, prop_schema: dict):

    try:
        if prop_schema.get('type') == 'integer':
            converted_value = int(value)
        elif prop_schema.get('type') == 'number':
            converted_value = float(value)
        elif prop_schema.get('type') == 'array':
            converted_value = ast.literal_eval(value)
        else:
            converted_value = value
    except ValueError:
        converted_value = value  # Let validation catch the error

    return converted_value


def validate_tsv_with_schema(tsv_file_path: str, schema_type: ValidationSchema):
    """
    Validates a TSV file (optionally gzipped) against a given JSON schema.
    Each row is validated as a JSON object (dict).
    Raises jsonschema.ValidationError if any row is invalid.
    """
    schema_path = files('mpralib.schemas').joinpath(schemaFilemap.get(schema_type))
    with schema_path.open('r', encoding='utf-8') as f:
        schema = json.load(f)

    open_func = gzip.open if is_compressed_file(tsv_file_path) else open
    mode = 'rt' if is_compressed_file(tsv_file_path) else 'r'

    with open_func(tsv_file_path, mode, encoding='utf-8') as tsvfile:
        reader = list(csv.DictReader(tsvfile, delimiter='\t'))
        for i, row in enumerate(tqdm.tqdm(reader, desc="Validating rows"), start=1):
            # Convert types according to schema
            for prop, prop_schema in schema.get('properties', {}).items():
                if prop in row and row[prop] != '':
                    if "anyOf" in prop_schema:
                        for anyOfProp_schema in prop_schema["anyOf"]:
                            row[prop] = _convert_row_value(row[prop], anyOfProp_schema)
                    else:
                        row[prop] = _convert_row_value(row[prop], prop_schema)
            try:
                jsonschema.validate(instance=row, schema=schema)
            except jsonschema.ValidationError as e:
                print(f"Row {i} invalid: {e.message}")
                raise e
            except Exception as e:
                print(f"Row {i} error: {e}")
                raise e
