module.exports = function (api) {
  api.cache(true);

  const plugins = ["@babel/plugin-proposal-class-properties"];
  const presets =  [
    [
      "@babel/preset-env",
      {
        "useBuiltIns": "entry"
      }
    ]
  ];

  return {
    presets,
    plugins
  };
}
