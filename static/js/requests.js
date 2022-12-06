var http = new XMLHttpRequest();

function getResponseFromURL(queryUrl) {
  http.open('GET', queryUrl, false);
  http.send(null);
  if (http.status === 200) {
    return http.responseText;
  } else {
    return null;
  }
}