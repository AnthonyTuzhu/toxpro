
function plotHeatMap(root, data) {


    // set the dimensions and margins of the graph
    var margin = {top: 30, right: 30, bottom: 30, left: 250},
      width = 750 - margin.left - margin.right,
      height = 750 - margin.top - margin.bottom;

    // append the svg object to the body of the page

    var svg = root.append("svg")
              .attr("width", width + margin.left + margin.right)
              .attr("height", height + margin.top + margin.bottom)
              .append("g")
              .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    function onlyUnique(value, index, self) {
    return self.indexOf(value) === index;
    }

    // unique CIDS
    var cids = data.cids;
    var aids = data.aids;
    // Labels of row and columns

    // Build X scales and axis:
    var x = d3.scaleBand()
      .range([ 0, width ])
      .domain(aids)
      .padding(0.01);
    svg.append("g")
      .attr("transform", "translate(0," + height + ")")
      .call(d3.axisBottom(x))
        .selectAll("text")
        .style("text-anchor", "end")
        .attr("dx", "-.8em")
        .attr("dy", ".15em")
        .attr("transform", "rotate(-65)")
        .style("font-size",  0);


    // Build X scales and axis:
    var y = d3.scaleBand()
      .range([ height, 0 ])
      .domain(cids)
      .padding(0.01);
    svg.append("g")
      .call(d3.axisLeft(y))
      .style("font-size",  0);

    // Build color scale
    var myColor = d3.scaleOrdinal()
      .range(["blue", "white", "red"])
      .domain([-1,0,1]);

    var dataPlot = data.cids.map(function(e, i) {
          return [e, data.aids[i], data.outcomes[i]];
        });

    svg.selectAll()
      .data(dataPlot)
      .enter()
      .append("rect")
      .attr("x", function(d, i) { return x(d[1]) })
      .attr("y", function(d, i) { return y(d[0]) })
      .attr("width", x.bandwidth() )
      .attr("height", y.bandwidth() )
      .style("fill", function(d) { return myColor(d[2])} )

}