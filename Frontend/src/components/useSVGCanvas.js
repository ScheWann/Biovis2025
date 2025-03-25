import React, { useEffect, useState, useLayoutEffect } from "react";
import * as d3 from "d3";

export default function useSVGCanvas(d3Container) {
  //takes a ref to a container, makes an svg over it, and returns the svg selection, size ,and a tooltip
  const [windowHeight, windowWidth] = useWindowSize();
  const [height, setHeight] = useState(0);
  const [width, setWidth] = useState(0);

  //todo: figue out what type a div is
  const [svg, setSvg] = useState();
  const [tTip, setTTip] = useState();

  useEffect(() => {
    if (d3Container.current) {
      d3.select(d3Container.current).selectAll("svg").remove();

      var h = d3Container.current.clientHeight;
      var w = d3Container.current.clientWidth;

      var canvas = d3
        .select(d3Container.current)
        .append("svg")
        .attr("class", "frameEntryD3")
        .attr("width", w)
        .attr("height", h);

      if (d3.select("body").select(".tooltip").empty()) {
        d3.select("body")
          .append("div")
          .attr("class", "tooltip")
          .style("visibility", "hidden");
      }
      var tip = d3.select("body").select(".tooltip");

      setHeight(h);
      setWidth(w);
      setSvg(canvas);
      setTTip(tip);
    }
  }, [d3Container.current, windowWidth, windowHeight]);

  return [svg, height, width, tTip];
}

function useWindowSize() {
  const [size, setSize] = useState([0, 0]);
  useLayoutEffect(() => {
    function updateSize() {
      setSize([window.innerWidth, window.innerHeight]);
    }
    window.addEventListener("resize", updateSize);
    updateSize();
    return () => window.removeEventListener("resize", updateSize);
  }, []);
  return size;
}
