from scraper import getdocumentation, geticon, extract, scrapeInputs



data = scrapeInputs()

def printInputs(classname,indent="",prefix=""):
    if classname not in data.keys(): return ""
    html = ""
    for input in data[classname]["inputs"]:
        if input["type"] == "group": continue
        if (input["type"] == "select"):
            id = prefix + input["string"]
            html     += indent + """<div class="card mb-3" >\n"""
            html     += indent + """  <div class="card-header">\n"""
            html     += indent + """    <h3 class="card-title mt-3 mb-1">{}</h3>\n""".format(id + ".type")
            html     += indent + """    <ul class="nav nav-tabs card-header-tabs" data-group="{}" role="tablist">\n""".format(id + "-tab")
            for cls in input["classes"]:
                name = cls.split("::")[-1].lower()
                html += indent + """      <li class="nav-item">\n"""
                html += indent + """        <a class="nav-link" id="{}" data-bs-toggle="tab" href="#{}" role="tab" >{}</a> \n""".format(name, input["string"] + "." + name + "-div", name)
                html += indent + """      </li>\n"""
            html +=     indent + """    </ul>\n"""
            html +=     indent + """    <input type="hidden" id="{}" name="{}" value="">\n""".format(id+"-id",id+".type")
            html +=     indent + """  </div>\n"""
            html +=     indent + """  <div class="card-body">\n"""
            html +=     indent + """    <div class="tab-content">\n"""
            for cls in input["classes"]:
                name = cls.split("::")[-1].lower()
                html += indent + """      <div class="tab-pane" id="{}" role="tabpanel">""".format(input["string"] + "." + name + "-div")
                html += indent + "<p>"
                html += data[cls]["documentation"]
                html += indent + "</p>"
                html += printInputs(cls,indent+"      ",
                                    prefix = input["string"] + "." + name + ".")
                html += indent + """      </div>\n"""
            html +=     indent + """    </div>\n"""
            html +=     indent + """  </div>\n"""
            html +=     indent + """</div>\n"""
            #html     += indent + """<input type='hidden' name='{}' value='[not completed yet - please enter manually]' >\n""".format(id+".type")
        elif (input["type"] == "queryclass"):
            id = prefix + input["string"]
            html     += indent + """<div class="card mb-3" >\n"""
            html     += indent + """  <div class="card-header">\n"""
            html     += indent + """    <h3 class="card-title mt-3 mb-1">{}</h3>\n""".format(id)
            html +=     indent + """    </ul>\n"""
            html +=     indent + """  </div>\n"""
            html +=     indent + """  <div class="card-body">\n"""
            html +=     indent + """    <div class="tab-content">\n"""
            cls = input["class"]
            cls = cls.split("<")[0] # remove template args for now
            if not cls in data.keys(): # if it's not an actual classname then
                cls = "::".join(classname.split("::")[:-1])+"::"+cls
            print(cls)
            name = cls.split("::")[-1].lower()
            html += printInputs(cls,indent+"      ",
                                prefix = input["string"] + ".")
            html +=     indent + """    </div>\n"""
            html +=     indent + """  </div>\n"""
            html +=     indent + """</div>\n"""
        elif (input["type"] == "querysubclass"):
            cls = input["class"]
            cls = cls.split("<")[0] # remove template args for now
            if not cls in data.keys(): # if it's not an actual classname then
                cls = "::".join(classname.split("::")[:-1])+"::"+cls
            print(cls)
            html += printInputs(cls,indent=indent, prefix=prefix)
        elif "validate" in input["type"]:
            id = prefix + input["string"]
            html += '<div class="input-group mb-3">'
            html += '  <div class="input-group-prepend">'
            html += indent + "    <span class='input-group-text'> {} </span>\n".format(id)
            html += indent + "  </div>\n"
            html += indent + """  <select class="form-select" id="{}"  name="{}">""".format(id,id)
            html += indent + """  <option selected disabled value="">Choose...</option>"""
            print(input["possibles"].split(','))
            for option in input["possibles"].replace('"','').split(','):
                html += indent + """  <option>{}</option>""".format(option)
            html += indent + """  </select>"""
            html += indent + "</div>"
        else:
            id = prefix + input["string"]
            html += '<div class="input-group mb-3">'
            html += '  <div class="input-group-prepend">'
            html += indent + """<span type="button" class="input-group-text" data-bs-toggle="popover" data-bs-title="{}" data-bs-content="{}">{}</span>""".format(input["file"].split("/")[-1]+":"+str(input["line"])
                                                                                                                                                                  ,input["doc"],id)
            #html += indent + "    <a class='input-group-text' data-bs-toggle='popover' data-content='documentation'> {} </a>\n".format(id)
            html += indent + "  </div>\n"
            required = ""
            if 'required' in input["type"]: required="required"
            placeholder = ""
            if input["type"] in ["query_default","queryarr_default"]: placeholder=input["default"]
            html += indent + "  <input type='text' id='{}' name='{}' class='form-control' {} placeholder='{}'>\n".format(
                id,id,required,placeholder)
            html += indent + "</div>"
    return html


html = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Dynamic Form with Nested Sections</title>
    <!-- Bootstrap CSS -->
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-T3c6CoIi6uLrA9TneNEoa7RxnatzjcDSCmG1MXxSR1GAsXEV/Dwwykc2MPK8M2HN" crossorigin="anonymous">
</head>
<body>
<style>
        .code-container {
            position: relative;
            background-color: #f8f9fa;
            border: 1px solid #ddd;
            border-radius: 5px;
            padding: 10px;
            font-family: monospace;
            white-space: pre-wrap;
            overflow-x: auto;
        }
        .variable { color: #007bff; } /* Blue for variables */
        .string { color: #28a745; }   /* Green for strings */
        .comment { color: #6c757d; }  /* Grey for comments */
        .copy-button {
            position: absolute;
            top: 10px;
            right: 10px;
        }
</style>
<div class='container mt-5 mb-5'> 
<h2 class='mb-4'>Integrator::Flame inputs</h2><form  id='form'>

"""

html += printInputs("Integrator::Flame")

html += r"""
</form>
<form>
  <h2 class='mb-4'>Alamo input file</h2><form  id='form'>
  <div class="form-group">
    <label for="code-container">The following contains only set values and is updated dynamically.</label>
    <textarea class="form-control" id="code-container" rows="10"></textarea>
  </div>
</form>
</div>
</body>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js" integrity="sha384-C6RzsynM9kWDrMNeT87bh95OGNyZPhcTNXj1NW7RuBCsyN/o0jlpcV8Qyq46cDfL" crossorigin="anonymous"></script>
<script src="https://cdn.jsdelivr.net/npm/@popperjs/core@2.11.8/dist/umd/popper.min.js"></script>
<script>


// update textbox with input file
function updateResult() {
    const form = document.getElementById('form');
    const inputs = form.querySelectorAll('input, select');
    var concatenatedString = "alamo.program = flame\n"
    inputs.forEach((input) => {
      if(input.value.trim() != '')
          concatenatedString += input.name + " = " + input.value + "\n";
    });
    document.getElementById("code-container").value = concatenatedString;
}


// Add event listeners to all inputs within the form (except the result textbox)
document.querySelectorAll("input, select").forEach(input => {
    input.addEventListener("input", updateResult); // Trigger on any input change
});


// JavaScript for live validation
document.addEventListener('DOMContentLoaded', function () {


    // set hidden variables when tabs are selected
    document.body.addEventListener('shown.bs.tab', function (event) {
        // Get the tab group based on data-group attribute
        const tabGroup = event.target.closest('ul[data-group]');
        if (tabGroup) {
            const groupName = tabGroup.getAttribute('data-group');
            console.log(groupName);
            const activeTabId = event.target.getAttribute('id');
            console.log(activeTabId);
            console.log("this",groupName.replace("-tab",".type"))
            // Update the hidden input corresponding to this tab group
            //const hiddenInput = document.querySelector(`input[name="activeTab${groupName.charAt(0).toUpperCase() + groupName.slice(1)}"]`);
            const hiddenInput = document.getElementById(groupName.replace("-tab","-id"));
            console.log(hiddenInput);
            if (hiddenInput) {
                console.log(hiddenInput.value);
                hiddenInput.value = activeTabId;
                console.log(hiddenInput.name,"=",hiddenInput.value);
                updateResult();
            }
        }
    });


    const form = document.getElementById('form');
    const req_inputs = form.querySelectorAll('input[required]');
    const all_inputs = form.querySelectorAll('input');

    // Initialize all popovers
    all_inputs.forEach((input) => {
        new bootstrap.Popover(input);
    });


    // Validate fields on page load
    req_inputs.forEach((input) => {
        if (input.value.trim() === '') {
            input.classList.add('is-invalid');
        }
    });

    req_inputs.forEach((input) => {
        input.addEventListener('input', function () {
            if (!input.value) {
                input.classList.add('is-invalid');
            } 
            else if (input.value.trim() === '') {
                input.classList.add('is-invalid');
            } else {
                input.classList.remove('is-invalid');
            }
        });
    });

    all_inputs.forEach((input) => {
        input.addEventListener('input', function () {
            if (!input.value) {
                input.classList.remove('is-valid');
            } 
            else if (input.value.trim() === '') {
                input.classList.remove('is-valid');
            } else {
                input.classList.add('is-valid');
            }
        });
    });


    form.addEventListener('submit', function (event) {
        let isValid = true;

        inputs.forEach((input) => {
            if (input.value.trim() === '') {
                input.classList.add('is-invalid');
                isValid = false;
            } else {
                input.classList.remove('is-invalid');
            }
        });

        if (!isValid) {
            event.preventDefault();
            event.stopPropagation();
        }
    });


    // Add event listener for all tabs across all containers
    document.querySelectorAll('.tab-container').forEach(container => {
        const hiddenInput = container.querySelector('.selectedTab');
        const tabs = container.querySelectorAll('.nav-tabs .nav-link');

        tabs.forEach(tab => {
            tab.addEventListener('shown.bs.tab', function (event) {
                const activeTabId = event.target.id; // Get the ID of the active tab
                hiddenInput.value = activeTabId; // Update the hidden input value
                console.log(activeTabId);
            });
        });
    });
});


const popoverTriggerList = document.querySelectorAll('[data-bs-toggle="popover"]')
const popoverList = [...popoverTriggerList].map(popoverTriggerEl => new bootstrap.Popover(popoverTriggerEl))

</script>

"""

f = open("Builder.html",'w')
f.write(html)
f.close()
#print(html)
