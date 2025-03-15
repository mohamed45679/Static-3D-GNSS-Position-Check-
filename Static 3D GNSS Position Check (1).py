#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import dash_table
import math
import numpy as np
import pandas as pd
from pyproj import Transformer
import io
from xhtml2pdf import pisa

# ------------------- Utility Functions -------------------

def dms_to_decimal(degrees, minutes, seconds):
    """Convert DMS to decimal degrees."""
    if None in [degrees, minutes, seconds]:
        return None
    try:
        degrees = float(degrees)
        minutes = float(minutes)
        seconds = float(seconds)
        # Basic range checks
        if not (0 <= degrees <= 360 and 0 <= minutes < 60 and 0 <= seconds < 60):
            return None
        return degrees + minutes / 60 + seconds / 3600
    except ValueError:
        return None

def decimal_to_dms(decimal_degrees):
    """Convert decimal degrees to DMS string."""
    if decimal_degrees is None:
        return "N/A"
    degrees = int(decimal_degrees)
    minutes = int((decimal_degrees - degrees) * 60)
    seconds =round(((decimal_degrees - degrees) * 60 - minutes) * 60,2)
    return f"{degrees:03d}° {minutes:02d}' {seconds:05.2f}''"

def normalize_angle(angle):
    """Normalize angle to be within [-180, 180]."""
    while angle > 180:
        angle -= 360
    while angle < -180:
        angle += 360
    return angle

def normalize_0_360(angle):
    """Normalize angle to be within [0, 360)."""
    return angle % 360

def dms_input(id_prefix):
    """
    Helper function to build 3 input fields for DMS (Degrees, Minutes, Seconds).
    """
    return dbc.Row([
        dbc.Col(
            dbc.Input(
                id=f"{id_prefix}-degrees",
                placeholder="Degrees",
                type="number",
                min=0,
                max=360,
                step=1,
                required=True
            ),
            width=4
        ),
        dbc.Col(
            dbc.Input(
                id=f"{id_prefix}-minutes",
                placeholder="Minutes",
                type="number",
                min=0,
                max=59,
                step=1,
                required=True
            ),
            width=4
        ),
        dbc.Col(
            dbc.Input(
                id=f"{id_prefix}-seconds",
                placeholder="Seconds",
                type="number",
                min=0,
                max=59.99,
                step=0.01,
                required=True
            ),
            width=4
        ),
    ], className="mb-2")

def calculate_grid_values(e1, n1, e2, n2, hemisphere):
    """
    Calculate both:
      - Grid Bearing (decimal degrees, plus DMS)
      - Grid Convergence (decimal degrees, plus DMS)
    """
    if None in [e1, n1, e2, n2, hemisphere]:
        return None, "N/A", None, "N/A"

    try:
        # تقديري لكيفية حساب منطقة UTM (قد تحتاج لتعديل حسب الدولة)
        utm_zone = int((e1 // 1000000) % 60) + 1
        south = (hemisphere.lower() == "south")

        # Transformer from UTM -> Geographic (EPSG:4326)
        transformer_utm_to_geo = Transformer.from_crs(
            f"+proj=utm +zone={utm_zone} +{'south' if south else ''} +ellps=WGS84",
            "epsg:4326",
            always_xy=True
        )

        lon1, lat1 = transformer_utm_to_geo.transform(e1, n1)
        lon2, lat2 = transformer_utm_to_geo.transform(e2, n2)

        # Grid Bearing
        delta_e = e2 - e1
        delta_n = n2 - n1
        bearing_rad = math.atan2(delta_e, delta_n)
        grid_bearing_dd = (math.degrees(bearing_rad) + 360) % 360

        # Grid Convergence
        central_meridian = utm_zone * 6 - 183
        grid_convergence_rad = math.atan(
            math.tan(math.radians(lon1 - central_meridian)) * math.sin(math.radians(lat1))
        )
        grid_convergence_dd = math.degrees(grid_convergence_rad)

        return (
            grid_bearing_dd,
            decimal_to_dms(grid_bearing_dd),
            grid_convergence_dd,
            decimal_to_dms(grid_convergence_dd)
        )
    except Exception as e:
        return None, "N/A", None, "N/A"

# ------------------- Dash App & Layout -------------------

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "Static 3D GNSS Position Check"

# تحضير 19 صفًا فارغًا لجدول Gyrocompass مع الأعمدة الجديدة
default_gyro_data = []
for i in range(1, 21):  # 20 rows
    default_gyro_data.append({
        "utc_time": "",
        "direction_deg": None,
        "direction_min": None,
        "direction_sec": None,
        "obs_distance": None,
        "obs_v_dist": None,  # Observed Vertical Distance
        "height": None,      # Target Height (في حال الحاجة)
        "dgps_e": None,      # DGPS Easting
        "dgps_n": None,      # DGPS Northing
        "dgps_h": None       # DGPS Height
    })

app.layout = dbc.Container([
    # --- Header with Title, Device Name Input and Logo ---
    dbc.Row([
        dbc.Col([
            html.H1("Static 3D GNSS Position Check", className="text-center"),
            dbc.Input(id="device-name", placeholder="Enter Device Name", type="text", style={"margin-top": "10px"})
        ], width=8),
        dbc.Col([
            # مثال لشعار أو صورة
            html.Img(src="https://via.placeholder.com/150", style={"max-width": "100%"})
        ], width=4, className="text-center")
    ], className="my-4"),

    # --- Job Details ---
    dbc.Card([
        dbc.CardHeader(html.H5("Job Details")),
        dbc.CardBody([
            dbc.Row([
                dbc.Col([html.H6("Job Number"), dbc.Input(id="job-number", required=True)], width=3),
                dbc.Col([html.H6("Job Description"), dbc.Input(id="job-description", required=True)], width=3),
                dbc.Col([html.H6("Client"), dbc.Input(id="client", required=True)], width=3),
                dbc.Col([html.H6("Party Chief"), dbc.Input(id="party-chief", required=True)], width=3),
            ]),
            dbc.Row([
                dbc.Col([html.H6("Surveyor"), dbc.Input(id="surveyor", required=True)], width=3),
                dbc.Col([html.H6("Wharf"), dbc.Input(id="wharf", required=True)], width=3),
                dbc.Col([html.H6("Vessel"), dbc.Input(id="vessel", required=True)], width=3),
                dbc.Col([
                    html.H6("Date"),
                    dcc.DatePickerSingle(id="date", placeholder="Select a date", display_format='YYYY-MM-DD')
                ], width=3),
            ]),
            dbc.Row([
                dbc.Col([
                    html.H6("Hemisphere"),
                    dbc.Select(
                        id="hemisphere",
                        options=[
                            {"label": "North", "value": "North"},
                            {"label": "South", "value": "South"},
                        ],
                        placeholder="Select Hemisphere",
                        required=True
                    )
                ], width=3),
            ])
        ])
    ], className="mb-4"),

    # --- Instrument Station ---
    dbc.Card([
        dbc.CardHeader(html.H5("Instrument Station")),
        dbc.CardBody([
            dbc.Row([
                dbc.Col([html.H6("Instrument Station Name"), dbc.Input(id="instrument-station", required=True)], width=3),
                dbc.Col([html.H6("Easting (m)"), dbc.Input(id="easting", type="number", min=0, required=True)], width=3),
                dbc.Col([html.H6("Northing (m)"), dbc.Input(id="northing", type="number", min=0, required=True)], width=3),
                dbc.Col([html.H6("AHD Height (m)"), dbc.Input(id="AHD-height", type="number", required=True)], width=3),
            ]),
            # الحقل الجديد: Instrument Height
            dbc.Row([
                dbc.Col([html.H6("Instrument Height (m)"), dbc.Input(id="instrument-height", type="number", required=False)], width=3),
            ]),
        ])
    ], className="mb-4"),

    # --- Backsight Station ---
    dbc.Card([
        dbc.CardHeader(html.H5("Backsight Station")),
        dbc.CardBody([
            dbc.Row([
                dbc.Col([html.H6("Backsight Station Name"), dbc.Input(id="backsight-station", required=True)], width=3),
                dbc.Col([html.H6("Backsight Easting (m)"), dbc.Input(id="backsight-easting", type="number", min=0, required=True)], width=3),
                dbc.Col([html.H6("Backsight Northing (m)"), dbc.Input(id="backsight-northing", type="number", min=0, required=True)], width=3),
                dbc.Col([html.H6("Backsight AHD Height (m)"), dbc.Input(id="backsight-AHD-height", type="number", required=True)], width=3),
            ]),
        ])
    ], className="mb-4"),

    # --- Calculated Grid Bearing (DMS) ---
    dbc.Card([
        dbc.CardHeader(html.H5("Calculated Grid Bearing (DMS)")),
        dbc.CardBody([
            html.Div(id="grid-bearing-output", style={"fontWeight": "bold", "fontSize": "16px"})
        ])
    ], className="mb-4"),

    # --- Calculated Grid Convergence (DMS) ---
    dbc.Card([
        dbc.CardHeader(html.H5("Calculated Grid Convergence (DMS)")),
        dbc.CardBody([
            html.Div(id="grid-convergence-output", style={"fontWeight": "bold", "fontSize": "16px"})
        ])
    ], className="mb-4"),

    # --- Backsight Observation (DMS) ---
    dbc.Card([
        dbc.CardHeader(html.H5("Backsight Observation (DMS)")),
        dbc.CardBody([
            dms_input("backsight-observation")
        ])
    ], className="mb-4"),

    # --- Gyrocompass Observations: 19 Fixed Rows ---
    dbc.Card([
        dbc.CardHeader(html.H5("Gyrocompass Observations (20 Rows)")),
        dbc.CardBody([
            dash_table.DataTable(
                id='gyro-table',
                columns=[
                    {"name": "UTC Time (hh:mm:ss)",        "id": "utc_time",     "type": "text",    "editable": True},
                    {"name": "Direction Deg",              "id": "direction_deg","type": "numeric", "editable": True},
                    {"name": "Direction Min",              "id": "direction_min","type": "numeric", "editable": True},
                    {"name": "Direction Sec",              "id": "direction_sec","type": "numeric", "editable": True},
                    {"name": "Observed Distance (m)",      "id": "obs_distance", "type": "numeric", "editable": True},
                    {"name": "Observed Vertical Distance (m)", "id": "obs_v_dist","type": "numeric", "editable": True},
                    {"name": "Height (m)",                 "id": "height",       "type": "numeric", "editable": True},
                    {"name": "DGPS Easting (m)",           "id": "dgps_e",       "type": "numeric", "editable": True},
                    {"name": "DGPS Northing (m)",          "id": "dgps_n",       "type": "numeric", "editable": True},
                    {"name": "DGPS Height (m)",            "id": "dgps_h",       "type": "numeric", "editable": True},
                ],
                data=default_gyro_data,
                editable=True,
                row_deletable=False,
                style_cell={
                    'minWidth': '100px', 'width': '150px', 'maxWidth': '200px',
                    'whiteSpace': 'normal',
                    'textAlign': 'left',
                },
                style_header={
                    'backgroundColor': 'rgb(230, 230, 230)',
                    'fontWeight': 'bold'
                },
            ),
        ])
    ], className="mb-4"),

    # --- Submit Button for Main Calculations ---
    dbc.Row([
        dbc.Col([
            dbc.Button("submit_Horizontal", id="submit-button", color="primary", className="mt-4", n_clicks=0)
        ], width=12, className="text-center"),
    ], className="mb-4"),

    # --- Display: Job Details, Results, and Statistical Analysis ---
    dbc.Row([
        dbc.Col([
            html.H4("Job Details"),
            html.Div(id="submitted-job-details")
        ], width=12),
    ], className="mb-4"),

    dbc.Row([
        dbc.Col([
            html.H4("Calculated Results"),
            html.Div(id="results-table")
        ], width=12),
    ], className="mb-4"),

    dbc.Row([
        dbc.Col([
            html.H4("Statistical Analysis"),
            html.Div(id="statistical-analysis", style={"fontWeight": "bold", "fontSize": "18px"})
        ], width=12),
    ], className="mb-4"),

    # --- Submit Button for Vertical Calculations ---
    dbc.Row([
        dbc.Col([
            dbc.Button("Submit Vertical Data", id="vertical-submit-button", color="info", className="mt-4", n_clicks=0)
        ], width=12, className="text-center"),
    ], className="mb-4"),

    # --- Vertical Calculations Results ---
    dbc.Row([
        dbc.Col([
            html.H4("Vertical Height Calculations"),
            html.Div(id="vertical-results")
        ], width=12),
    ], className="mb-4"),

    # --- Export to PDF Button & Download Component ---
    dbc.Row([
        dbc.Col([
            dbc.Button("Export to PDF", id="export-pdf-button", color="secondary", className="mt-4"),
            dcc.Download(id="download-pdf")
        ], width=12, className="text-center")
    ], className="mb-4"),

], fluid=True)

# ------------------- Callbacks -------------------

# 1) تحديث وعرض قيم Grid Bearing و Grid Convergence
@app.callback(
    [Output("grid-bearing-output", "children"),
     Output("grid-convergence-output", "children")],
    [
        Input("easting", "value"),
        Input("northing", "value"),
        Input("backsight-easting", "value"),
        Input("backsight-northing", "value"),
        Input("hemisphere", "value")
    ]
)
def update_grid_and_convergence(e1, n1, e2, n2, hemisphere):
    gb_dd, gb_dms, gc_dd, gc_dms = calculate_grid_values(e1, n1, e2, n2, hemisphere)
    return gb_dms, gc_dms

# 2) المعالجة الرئيسية للبيانات وعرض النتائج
@app.callback(
    [
        Output("results-table", "children"),
        Output("statistical-analysis", "children"),
        Output("submitted-job-details", "children"),
    ],
    Input("submit-button", "n_clicks"),
    [
        State("job-number", "value"),
        State("job-description", "value"),
        State("client", "value"),
        State("party-chief", "value"),
        State("surveyor", "value"),
        State("wharf", "value"),
        State("vessel", "value"),
        State("date", "date"),
        State("hemisphere", "value"),
        State("instrument-station", "value"),
        State("easting", "value"),
        State("northing", "value"),
        State("AHD-height", "value"),
        State("instrument-height", "value"),  # الحقل الجديد
        State("backsight-station", "value"),
        State("backsight-easting", "value"),
        State("backsight-northing", "value"),
        State("backsight-AHD-height", "value"),
        # Backsight Observation DMS
        State("backsight-observation-degrees", "value"),
        State("backsight-observation-minutes", "value"),
        State("backsight-observation-seconds", "value"),
        # Table
        State('gyro-table', 'data'),
        State('gyro-table', 'columns'),
        # We also need the displayed Grid Bearing & Convergence
        State("grid-bearing-output", "children"),
        State("grid-convergence-output", "children"),
        # Device name
        State("device-name", "value")
    ]
)
def submit_data(
    n_clicks,
    job_number, job_description, client, party_chief, surveyor,
    wharf, vessel, date, hemisphere,
    instrument_station, e1, n1, ahd1,
    instrument_height,
    backsight_station, e2, n2, ahd2,
    bs_obs_deg, bs_obs_min, bs_obs_sec,
    gyro_data, gyro_columns,
    grid_bearing_str, grid_convergence_str,
    device_name
):
    if n_clicks == 0:
        return "", "", ""

    # التحقق الأساسي من المدخلات المطلوبة
    required_inputs = [
        job_number, job_description, client, party_chief, surveyor,
        wharf, vessel, date, hemisphere,
        instrument_station, e1, n1, ahd1,
        backsight_station, e2, n2, ahd2,
        bs_obs_deg, bs_obs_min, bs_obs_sec
    ]
    if any(val is None for val in required_inputs):
        return html.Div("Please fill all required fields."), "", ""

    # تحويل الـ Backsight Observation من DMS إلى درجة عشرية
    backsight_obs_dd = dms_to_decimal(bs_obs_deg, bs_obs_min, bs_obs_sec)
    if backsight_obs_dd is None:
        return html.Div("Backsight Observation (DMS) is invalid."), "", ""

    # إعادة حساب Grid Bearing & Convergence كقيم عشرية
    gb_dd, gb_dms, gc_dd, gc_dms = calculate_grid_values(e1, n1, e2, n2, hemisphere)
    if gb_dd is None or gc_dd is None:
        return html.Div("Could not compute Grid Bearing / Convergence."), "", ""

    # سنجمع النتائج في قائمة
    results_rows = []
    miscloses = []
    total_rows = 0
    successful_rows = 0
    failed_rows = 0

    for row in gyro_data:
        time_val = (row.get('utc_time') or "").strip()
        direction_deg = row.get('direction_deg')
        direction_min = row.get('direction_min')
        direction_sec = row.get('direction_sec')
        obs_distance = row.get('obs_distance')
        dgps_e = row.get('dgps_e')
        dgps_n = row.get('dgps_n')
        total_rows += 1
        # إذا كان السطر فارغًا تمامًا نتجاهله
        if not time_val and not direction_deg and not direction_min and not direction_sec and not obs_distance:
            continue

        # التحقق من صيغة الوقت
        try:
            hh, mm, ss = map(int, time_val.split(':'))
            if not (0 <= hh < 24 and 0 <= mm < 60 and 0 <= ss < 60):
                raise ValueError
        except:
            failed_rows += 1
            continue

        # تحويل DMS إلى درجة عشرية
        obs_dir_dd = dms_to_decimal(direction_deg, direction_min, direction_sec)
        if obs_dir_dd is None or obs_distance is None or obs_distance <= 0:
            failed_rows += 1
            continue

        # حساب plane bearing بالدرجات
        plane_bearing_dd = gb_dd + (obs_dir_dd - backsight_obs_dd)
        plane_bearing_dd = normalize_0_360(plane_bearing_dd)
        plane_bearing_dms = decimal_to_dms(plane_bearing_dd)

        # حساب الإحداثيات الناتجة عن هذه المسافة والزاوية
        plane_bearing_rad = math.radians(plane_bearing_dd)
        calc_e = e1 + obs_distance * math.sin(plane_bearing_rad)
        calc_n = n1 + obs_distance * math.cos(plane_bearing_rad)

        # حساب Linear Misclose إذا وجد DGPS Easting و DGPS Northing
        linear_misclose = None
        if dgps_e is not None and dgps_n is not None:
            try:
                diff_e = calc_e - dgps_e
                diff_n = calc_n - dgps_n
                linear_misclose = math.sqrt(diff_e**2 + diff_n**2)
                miscloses.append(linear_misclose)
            except:
                linear_misclose = None

        results_rows.append({
            "UTC Time": time_val,
            "Plane Bearing (DMS)": plane_bearing_dms,
            "Obs Dist (m)": f"{obs_distance:.3f}",
            "Calc Easting": f"{calc_e:.3f}",
            "Calc Northing": f"{calc_n:.3f}",
            "DGPS Easting": f"{dgps_e:.3f}" if dgps_e is not None else "N/A",
            "DGPS Northing": f"{dgps_n:.3f}" if dgps_n is not None else "N/A",
            "Linear Misclose (m)": f"{linear_misclose:.3f}" if linear_misclose is not None else "N/A"
        })
        successful_rows += 1

    # حساب المتوسط والانحراف المعياري للـ Misclose
    mean_misclose = None
    std_misclose = None
    if miscloses:
        mean_misclose = np.mean(miscloses)
        std_misclose = np.std(miscloses)

    # بناء جدول تفاصيل المشروع
    job_details = {
        "Job Number": job_number,
        "Job Description": job_description,
        "Client": client,
        "Party Chief": party_chief,
        "Surveyor": surveyor,
        "Wharf": wharf,
        "Vessel": vessel,
        "Date": date,
        "Hemisphere": hemisphere,
        "Device Name": device_name if device_name else "N/A",
        "Instrument Station": instrument_station,
        "Instrument Easting": e1,
        "Instrument Northing": n1,
        "Instrument AHD": ahd1,
        "Instrument Height": instrument_height if instrument_height is not None else "N/A",
        "Backsight Station": backsight_station,
        "BS Easting": e2,
        "BS Northing": n2,
        "BS AHD": ahd2,
        "Calculated Grid Bearing (DMS)": gb_dms or "N/A",
        "Calculated Grid Convergence (DMS)": gc_dms or "N/A"
    }
    job_details_table = html.Table([
        html.Tbody([
            html.Tr([html.Th(k), html.Td(str(v))]) for k, v in job_details.items()
        ])
    ], className="table table-bordered table-striped")

    # بناء جدول النتائج
    if results_rows:
        columns = [
            "UTC Time",
            "Plane Bearing (DMS)",
            "Obs Dist (m)",
            "Calc Easting",
            "Calc Northing",
            "DGPS Easting",
            "DGPS Northing",
            "Linear Misclose (m)"
        ]
        thead = html.Thead([
            html.Tr([html.Th(col) for col in columns])
        ])
        tbody = html.Tbody([
            html.Tr([
                html.Td(row[col]) for col in columns
            ]) for row in results_rows
        ])
        results_table = html.Table([thead, tbody], className="table table-bordered table-striped")
    else:
        results_table = html.Div("No valid Gyrocompass Observations to display.")

    # التحليل الإحصائي
    stats = [
        f"Total Observations: {total_rows}",
        f"Successful Calculations: {successful_rows}",
        f"Failed Calculations: {failed_rows}"
    ]
    if mean_misclose is not None and std_misclose is not None:
        stats.append(f"Mean(Linear Misclose): {mean_misclose:.2f} m")
        stats.append(f"Std Dev(Linear Misclose): {std_misclose:.2f} m")

    stat_div = html.Div([html.Ul([html.Li(s) for s in stats])])

    return results_table, stat_div, job_details_table

# 3) Callback لحفظ التقرير بصيغة PDF
def convert_html_to_pdf(html_content):
    """Convert raw HTML to PDF bytes using xhtml2pdf (pisa)."""
    pdf_file = io.BytesIO()
    pisa_status = pisa.CreatePDF(io.StringIO(html_content), dest=pdf_file)
    if pisa_status.err:
        return None
    pdf_file.seek(0)
    return pdf_file.getvalue()

@app.callback(
    Output("download-pdf", "data"),
    Input("export-pdf-button", "n_clicks"),
    State("job-number", "value"),
    State("job-description", "value"),
    State("client", "value"),
    State("party-chief", "value"),
    State("surveyor", "value"),
    State("wharf", "value"),
    State("vessel", "value"),
    State("date", "date"),
    State("hemisphere", "value"),
    State("instrument-station", "value"),
    State("easting", "value"),
    State("northing", "value"),
    State("AHD-height", "value"),
    State("instrument-height", "value"),
    State("backsight-station", "value"),
    State("backsight-easting", "value"),
    State("backsight-northing", "value"),
    State("backsight-AHD-height", "value"),
    State("grid-bearing-output", "children"),
    State("grid-convergence-output", "children"),
    State("results-table", "children"),
    State("statistical-analysis", "children"),
    prevent_initial_call=True
)
def export_to_pdf(n_clicks,
                  job_number, job_description, client, party_chief, surveyor,
                  wharf, vessel, date, hemisphere,
                  instrument_station, easting, northing, AHD_height,
                  instrument_height,
                  backsight_station, backsight_easting, backsight_northing, backsight_AHD_height,
                  grid_bearing_output, grid_convergence_output,
                  results_table, stat_analysis):

    if not n_clicks:
        return None

    # بناء محتوى الـ HTML للتقرير
    html_content = f"""
    <html>
      <head>
        <meta charset="utf-8">
        <title>Survey Report</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 20px;
                line-height: 1.6;
                font-size: 12px;
            }}
            h1, h2 {{
                text-align: center;
                color: #333;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
                page-break-inside: auto;
            }}
            th, td {{
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
                word-wrap: break-word;
            }}
            th {{
                background-color: #f2f2f2;
            }}
        </style>
      </head>
      <body>
        <h1>Survey Report</h1>

        <h2>Job Details</h2>
        <table>
          <tr><th>Job Number</th><td>{job_number}</td></tr>
          <tr><th>Job Description</th><td>{job_description}</td></tr>
          <tr><th>Client</th><td>{client}</td></tr>
          <tr><th>Party Chief</th><td>{party_chief}</td></tr>
          <tr><th>Surveyor</th><td>{surveyor}</td></tr>
          <tr><th>Wharf</th><td>{wharf}</td></tr>
          <tr><th>Vessel</th><td>{vessel}</td></tr>
          <tr><th>Date</th><td>{date}</td></tr>
          <tr><th>Hemisphere</th><td>{hemisphere}</td></tr>
          <tr><th>Instrument Station</th><td>{instrument_station}</td></tr>
          <tr><th>Instrument E</th><td>{easting}</td></tr>
          <tr><th>Instrument N</th><td>{northing}</td></tr>
          <tr><th>Instrument AHD</th><td>{AHD_height}</td></tr>
          <tr><th>Instrument Height</th><td>{instrument_height if instrument_height else "N/A"}</td></tr>
          <tr><th>Backsight Station</th><td>{backsight_station}</td></tr>
          <tr><th>Backsight E</th><td>{backsight_easting}</td></tr>
          <tr><th>Backsight N</th><td>{backsight_northing}</td></tr>
          <tr><th>Backsight AHD</th><td>{backsight_AHD_height}</td></tr>
          <tr><th>Calculated Grid Bearing (DMS)</th><td>{grid_bearing_output}</td></tr>
          <tr><th>Calculated Grid Convergence (DMS)</th><td>{grid_convergence_output}</td></tr>
        </table>

        <h2>Calculated Results</h2>
        <div>{results_table if results_table else "<p>No Results Available</p>"}</div>

        <h2>Statistical Analysis</h2>
        <div>{stat_analysis if stat_analysis else "<p>No Statistical Analysis</p>"}</div>
      </body>
    </html>
    """

    pdf_bytes = convert_html_to_pdf(html_content)
    if pdf_bytes:
        return dcc.send_bytes(pdf_bytes, "survey_report.pdf")
    return None

# 4) Callback للمعالجة العمودية (Vertical Calculations)
@app.callback(
    Output("vertical-results", "children"),
    Input("vertical-submit-button", "n_clicks"),
    [
        State("AHD-height", "value"),
        State("instrument-height", "value"),
        State("gyro-table", "data")
    ]
)
def calculate_vertical_data(n_clicks, ahd_height, instrument_height, gyro_data):
    if n_clicks == 0:
        return ""
    if ahd_height is None or instrument_height is None:
        return html.Div("Please provide both AHD Height and Instrument Height for vertical calculations.")

    results_rows = []
    vertical_miscloses = []
    for row in gyro_data:
        time_val = (row.get("utc_time") or "").strip()
        # الحصول على القيمة المقاسة للمسافة العمودية (Observed Vertical Distance)
        try:
            obs_v = float(row.get("obs_v_dist"))
            obs_v2 = float(row.get("height"))
        except (TypeError, ValueError):
            continue
        # حساب calculated coordinate height وفقًا للصيغة الجديدة:
        # Calculated Coord. Height = Benchmark Height + Instrument Height + Observed Vertical Distance
        calc_height = ahd_height + instrument_height + obs_v+ obs_v2

        # محاولة حساب Vertical Misclose إذا كانت قيمة DGPS Height متوفرة
        dgps_h = row.get("dgps_h")
        vertical_misclose = None
        if dgps_h is not None:
            try:
                dgps_h = float(dgps_h)
                vertical_misclose = abs(calc_height - dgps_h)
                vertical_miscloses.append(vertical_misclose)
            except (TypeError, ValueError):
                vertical_misclose = None

        results_rows.append({
            "UTC Time": time_val,
            "Benchmark Height": f"{ahd_height:.3f}",
            "Instrument Height": f"{instrument_height:.3f}",
            "Observed Vertical Dist.": f"{obs_v:.3f}",
            "Calculated Co_ord. Height": f"{calc_height:.3f}",
            "DGPS Height": f"{float(dgps_h):.3f}" if dgps_h is not None and vertical_misclose is not None else "N/A",
            "Vertical Misclose": f"{vertical_misclose:.3f}" if vertical_misclose is not None else "N/A"
        })

    # بناء جدول النتائج العمودية
    if results_rows:
        columns = ["UTC Time", "Benchmark Height", "Instrument Height", "Observed Vertical Dist.", "Calculated Co_ord. Height", "DGPS Height", "Vertical Misclose"]
        table_header = html.Thead(html.Tr([html.Th(col) for col in columns]))
        table_body = html.Tbody([
            html.Tr([html.Td(row[col]) for col in columns]) for row in results_rows
        ])
        table = html.Table([table_header, table_body], className="table table-bordered table-striped")
    else:
        table = html.Div("No valid vertical data to display.")

    # حساب المتوسط والانحراف المعياري للـ Vertical Misclose
    stats = []
    if vertical_miscloses:
        mean_vm = np.mean(vertical_miscloses)
        std_vm = np.std(vertical_miscloses)
        stats = [
            html.Li(f"Mean Vertical Misclose: {mean_vm:.3f} m"),
            html.Li(f"Std Dev Vertical Misclose: {std_vm:.3f} m")
        ]
    stats_div = html.Div([html.Ul(stats)])
    return html.Div([table, stats_div])

# ------------------- Run Server -------------------
if __name__ == "__main__":
    app.run_server(debug=True)

