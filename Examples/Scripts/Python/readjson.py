import json

def read_json(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            data = json.load(file)
            return data
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except json.JSONDecodeError:
        print(f"Error: The file '{file_path}' is not a valid JSON file.")

def save_json(file_path, data):
    try:
        with open(file_path, 'w', encoding='utf-8') as file:
            json.dump(data, file, indent=4)
    except Exception as e:
        print(f"Error: Unable to save JSON file. {e}")

# Example usage
if __name__ == "__main__":
    file_path = "geometry-map.json"  # Change this to your JSON file path
    json_data = read_json(file_path)
    if json_data:
        surface = json_data["Surfaces"]["entries"]
        for auto in surface:
            if "approach" in auto:
                if auto["approach"] == 2:
                    auto["value"]["material"]["binUtility"]["binningdata"][1]["bins"] = 100
                    auto["value"]["material"]["mapMaterial"]=True
                elif auto["layer"] == 20 and auto["approach"] == 1:
                    auto["value"]["material"]["binUtility"]["binningdata"][1]["bins"] = 100
                    auto["value"]["material"]["mapMaterial"]=True
    save_json(file_path, json_data)
