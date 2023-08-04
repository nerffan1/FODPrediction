#!/usr/bin/python3

import requests

resp = requests.get("https://www.nist.gov/sites/default/files/data.json")
data = resp.json()
i = 0
print(data["@type"])
print(data["describedBy"])