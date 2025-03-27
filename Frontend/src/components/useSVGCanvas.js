import React, { useEffect, useState, useLayoutEffect } from "react";
import * as d3 from "d3";

export default function useSVGCanvas(d3Container) {
  const [windowHeight, windowWidth] = useWindowSize();
  const [height, setHeight] = useState(0);
  const [width, setWidth] = useState(0);
  const [svg, setSvg] = useState();

  useEffect(() => {
    if (d3Container.current) {
      // Clear any existing SVG
      d3.select(d3Container.current).selectAll("svg").remove();

      // Get container dimensions
      const h = d3Container.current.clientHeight;
      const w = d3Container.current.clientWidth;

      // Create new SVG canvas
      const canvas = d3
        .select(d3Container.current)
        .append("svg")
        .attr("class", "d3-canvas")
        .attr("width", w)
        .attr("height", h);

      // Update state
      setHeight(h);
      setWidth(w);
      setSvg(canvas);
    }
  }, [d3Container.current, windowWidth, windowHeight]);

  return [svg, height, width];
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
