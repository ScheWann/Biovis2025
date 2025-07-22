import { LineChart } from "./LineChart";

export const GeneExpression = () => {
  return (
    <div
      style={{
        display: "flex",
        flexDirection: "column",
      }}
    >
      <LineChart
        data={[2, 7, 4, 5]}
        xLabel="Index"
        yLabel="Value"
        title="Numbers"
      />
      <LineChart
        data={[
          { name: "A", val: 10 },
          { name: "B", val: 20 },
        ]}
        xAccessor={(d) => d.name}
        yAccessor={(d) => d.val}
        xLabel="Category"
        yLabel="Score"
        title="Objects"
      />
    </div>
  );
};
