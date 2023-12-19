import time
from datetime import datetime
import Fluigent.SDK as fluigent
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Update the graph with new data as it is collected
def update_plot(frame, data_time, data_pressure, line, ax, csv):
    # data_time[0] is the actual start value. The first data added to the csv is data_time[1]
    current_time = time.perf_counter() - data_time[0]
    measurement = fluigent.fgt_get_sensorValue(0)
    print(f'{measurement:.4f}')

    data_time.append(current_time)
    data_pressure.append(measurement)

    csv.writelines(f'{current_time:.4f},{measurement:.4f}\n')

    # Ignore the first value in the data (raw time.perf_counter() output. i.e. =/= 0)
    line.set_data(data_time[1:], data_pressure[1:])

    # Change the dimensions of the graph to fit the newly added data
    ax.relim()
    ax.autoscale_view()

    # Scroll graph showing 600 data at once
    max_data = 600
    if len(data_time) < max_data:
        ax.set_xlim(0, data_time[-1])
    else:
        ax.set_xlim(data_time[- max_data], data_time[-1])

    # # Show all data
    # ax.set_xlim(0, data_time[-1])

def main():
    data_time, data_pressure = [], []

    # # Create a simulated IPS instrument for dev purposes
    # fluigent.fgt_create_simulated_instr(5, 12345, 0, [11])

    # Detect any connected fluigent devices
    # Make sure they are not plugged in via a usb splitter. This causes initialisation issues.
    serial_numbers, types = fluigent.fgt_detect()

    for i, serial_num in enumerate(serial_numbers):
        print(f'Detected instrument at index: {i}, ControllerSN: {serial_num}, type: {str(types[i])}')

    # Initialise the devices in serial_numbers
    fluigent.fgt_init(serial_numbers)
    controllerInfoArray = fluigent.fgt_get_controllersInfo()

    for i, controllerInfo in enumerate(controllerInfoArray):
        print(f'Controller info at index: {i}')
        print(controllerInfo)

        # Calibrate readings to 0 using water (i) as a reference
        # Other solvents possible. See fgt_SENSOR_CALIBRATION in 'Fluigent SDK.pdf' Chapter 5
        fluigent.fgt_set_sensorCalibration(i, 1)

        title_string = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        with open(f'{title_string}.csv', 'w') as output_csv:
            current_time = time.perf_counter()
            measurement = fluigent.fgt_get_sensorValue(i)
            unit = fluigent.fgt_get_sensorUnit(i)
            output_csv.writelines(f'Time (s),Pressure ({unit})\n')

            fig, ax = plt.subplots()
            line, = ax.plot([], [], label='Pressure')
            ax.set_ylim(-10, 2232.0811)

            plt.title(title_string.replace('_', ' '))
            plt.xlabel('Time (s)')
            plt.ylabel(f'Pressure ({unit})')

            # Add raw time.perf_counter() value for first measurement
            data_time.append(current_time)
            data_pressure.append(measurement)

            # Write 0 values for time and pressure to signify the beginning of the outputted data
            output_csv.writelines(f'{0:.4f},{0:.4f}\n')
            data_time.append(0)
            data_pressure.append(0)

            # Update the data and graph every 500 ms 
            animation = FuncAnimation(fig,
                                      update_plot,
                                      fargs=(data_time,
                                             data_pressure,
                                             line,
                                             ax,
                                             output_csv),
                                      interval=500,
                                      cache_frame_data=False)
            
            # plt.legend()
            plt.show()

            fluigent.fgt_close()

if __name__ == '__main__':
    main()
