function addPredictorToolButton() {

            var dropdownDiv = $("#dropdown-page-tools");

            var button = $("                      <button class=\"btn btn-info\" style=\"margin-left: 10px\" type=\"button\" id=\"dropdownMenuButton\" data-toggle=\"dropdown\" aria-haspopup=\"true\" aria-expanded=\"false\">\n" +
                "                        <i class=\"fas fa-chevron-circle-down\"></i>\n" +
                "                        <span style=\"margin-left: 5px\">Actions</span>\n" +
                "                      </button>");

            var dropdown = $("<div class=\"dropdown-menu\" aria-labelledby=\"dropdownMenuButton\"></div>")


            var downloadPRofileLink = $("<a onclick=\"downloadPredictions()\" class=\"dropdown-item\"" +
                " href=\"#\"><i class=\"fas fa-download\"></i><span style=\"margin-left: 10px\">Download predictions</span></a>");

            dropdown.append(downloadPRofileLink);

            dropdownDiv.append(button);
            dropdownDiv.append(dropdown);
}